"""
索引构建器 (Incremental Version)
实现 Incremental DAG 构建逻辑，支持 Multi-Parent
解决 O(N^2) 内存问题，支持流式插入
"""

import uuid
from typing import List, Dict, Union, Set
from .structure import WorldNode, BioSequence, GenomePointer, R_SW, R_MW, R_LW
from .tools import compute_distance

class BioGeometryIndexBuilder:
    """索引构建器主类 (增量插入版)"""
    
    def __init__(self):
        """初始化索引构建器"""
        # 层级存储：Layer 3 (LW) 是入口
        self.layers: Dict[int, List[WorldNode]] = {
            1: [],  # Layer 1: SW
            2: [],  # Layer 2: MW
            3: []   # Layer 3: LW (Root Candidates)
        }
        
        # 半径配置
        self.radius_config = {
            1: R_SW,
            2: R_MW,
            3: R_LW
        }
        
        # 统计信息
        self.stats = {
            'added_sequences': 0,
            'created_nodes': {1: 0, 2: 0, 3: 0}
        }
    
    def build(self, raw_sequences: List[BioSequence]):
        """
        构建索引的主入口 (流式处理)
        
        Args:
            raw_sequences: 原始序列列表 (可以是生成器)
        """
        print(f"[Build] Starting Incremental Build for {len(raw_sequences)} sequences...")
        
        for i, seq in enumerate(raw_sequences):
            self.add_sequence(seq)
            
            if (i + 1) % 100 == 0:
                print(f"  Processed {i + 1} sequences...", end='\r')
        
        print(f"\n[Build] Completed. Added {self.stats['added_sequences']} sequences.")
        self._print_summary()

    def add_sequence(self, new_sequence: BioSequence):
        """
        插入单个新序列 (伪代码: InsertSequence)
        """
        self.stats['added_sequences'] += 1
        
        # 1. 寻找能包裹该序列的现有 SW 节点
        # SearchForParents(new_sequence, IndexRoot, Layer=SW)
        # 对于 Sequence，查询半径为 0
        valid_sws = self._search_candidates(
            query_center=new_sequence, 
            query_radius=0, 
            target_layer=1
        )

        if valid_sws:
            # Case 1: 找到了一个或多个能容纳它的 SW
            # 执行 DAG 多重插入
            for sw_node in valid_sws:
                sw_node.children.append(new_sequence)
                sw_node.data_count += 1
            return
        else:
            # Case 2: 没有任何旧 SW 能包含它
            # 创建新的 SW 节点
            new_sw_node = WorldNode(
                center_ptr=new_sequence, # 这里直接用 sequence，实际工程中可能是 Pointer
                radius=R_SW,
                layer_level=1
            )
            # 把自己加进去作为第一个孩子
            new_sw_node.children.append(new_sequence)
            new_sw_node.data_count = 1
            
            # 注册到 Layer 1 列表 (用于统计和调试)
            self.layers[1].append(new_sw_node)
            self.stats['created_nodes'][1] += 1
            
            # 递归向上传递 (InsertNodeUpwards)
            self._insert_node_upwards(new_sw_node)

    def _insert_node_upwards(self, child_node: WorldNode):
        """
        递归向上插入节点 (伪代码: InsertNodeUpwards)
        """
        current_layer = child_node.layer
        
        # 如果已经是顶层 (LW)，直接加入根列表并结束
        if current_layer == 3:
            self.layers[3].append(child_node)
            self.stats['created_nodes'][3] += 1
            return

        target_parent_layer = current_layer + 1
        parent_radius = self.radius_config[target_parent_layer]

        # 1. 在上一层寻找能包裹 child_node 的父节点
        # 判定标准: 严格包含 (Strict Inclusion)
        # Dist(Parent, Child) + Child.R <= Parent.R
        valid_parents = self._search_candidates(
            query_center=child_node.center_ptr,
            query_radius=child_node.radius,
            target_layer=target_parent_layer
        )

        if valid_parents:
            # Case A: 找到了父节点 (DAG 连接)
            for parent in valid_parents:
                parent.children.append(child_node)
            return
        else:
            # Case B: 找不到父节点
            # 创建新的父节点，以子节点的中心为中心
            new_parent_node = WorldNode(
                center_ptr=child_node.center_ptr,
                radius=parent_radius,
                layer_level=target_parent_layer
            )
            new_parent_node.children.append(child_node)
            
            # 注册到对应层级列表
            self.layers[target_parent_layer].append(new_parent_node)
            self.stats['created_nodes'][target_parent_layer] += 1
            
            # 继续递归
            self._insert_node_upwards(new_parent_node)

    def _search_candidates(self, query_center, query_radius: int, target_layer: int) -> List[WorldNode]:
        """
        自顶向下搜索符合几何约束的节点 (伪代码: SearchForParents)
        
        Args:
            query_center: 查询中心 (BioSequence 或 GenomePointer)
            query_radius: 查询对象的半径 (Sequence=0, SW=5, MW=15)
            target_layer: 目标层级 (1, 2, or 3)
            
        Returns:
            符合 Strict Inclusion 的目标层节点列表
        """
        # 从顶层 LW (Layer 3) 开始搜索
        # 如果目标就是 Layer 3，直接在 Layer 3 中搜
        # 如果目标是 Layer 1，则路径为 L3 -> L2 -> L1
        
        current_candidates = self.layers[3] # IndexRoot.LWs
        current_layer_level = 3
        
        # 1. 逐层下钻 (Drill Down)
        while current_layer_level > target_layer:
            next_candidates = []
            
            # 遍历当前层候选者
            for node in current_candidates:
                dist = compute_distance(node.center_ptr, query_center)
                
                # 粗筛 (Pruning): 三角不等式
                # 如果 Dist(Q, Node) > Node.R + Q.R，则两球分离，不可能包含
                # 只有重叠或相切才进去看
                if dist <= node.radius + query_radius:
                    # 只有 WorldNode 才有 children 列表指向下一层节点
                    # 注意：Layer 1 的 children 是 Sequence，不能 drill down，但循环条件 current > target 保证了不会对 L1 drill down
                    next_candidates.extend([c for c in node.children if isinstance(c, WorldNode)])
            
            # 去重 (DAG 中一个节点可能有多个父节点，防止重复加入 next_candidates)
            # WorldNode 需要有唯一 ID，这里简单用 set 转换
            unique_candidates = {n.node_id: n for n in next_candidates}
            current_candidates = list(unique_candidates.values())
            
            current_layer_level -= 1
            
            # 如果中间某一层被剪枝完了，直接返回空
            if not current_candidates:
                return []

        # 2. 最终检查 (Final Check in Target Layer)
        final_matches = []
        for node in current_candidates:
            dist = compute_distance(node.center_ptr, query_center)
            
            # 严格包含判定 (Strict Inclusion)
            # 父球必须完全包裹子球: Dist + Child.R <= Parent.R
            if dist + query_radius <= node.radius:
                final_matches.append(node)
                
        return final_matches

    def _print_summary(self):
        """打印构建统计"""
        print(f"  Layer 1 (SW): {len(self.layers[1])} nodes")
        print(f"  Layer 2 (MW): {len(self.layers[2])} nodes")
        print(f"  Layer 3 (LW): {len(self.layers[3])} nodes")
        
        if self.layers[1] and self.layers[2]:
            dag_links = sum(len(mw.children) for mw in self.layers[2])
            redundancy = (dag_links / len(self.layers[1]) - 1) * 100
            print(f"  DAG Redundancy: {redundancy:.2f}% (Avg parents per SW: {dag_links/len(self.layers[1]):.2f})")

    def get_statistics(self) -> Dict:
        """
        获取完整的索引统计信息
        """
        # 1. 基础计数
        raw_count = self.stats['added_sequences']  # O(1) 获取
        sw_count = len(self.layers[1])
        mw_count = len(self.layers[2])
        lw_count = len(self.layers[3])

        # 2. 压缩率计算 (防止除以零)
        comp_sw = raw_count / sw_count if sw_count > 0 else 0
        comp_mw = sw_count / mw_count if mw_count > 0 else 0
        comp_lw = mw_count / lw_count if lw_count > 0 else 0

        # 3. DAG 冗余度计算 (O(M), M=MW节点数，非常快)
        # 逻辑: 统计所有 MW 指向 SW 的指针总数 vs 唯一的 SW 节点数
        dag_redundancy = 0
        if sw_count > 0 and mw_count > 0:
            # 统计 Layer 2 (MW) 所有的子节点引用数
            # 注意: 在 Layer 2, child 就是 Layer 1 的 WorldNode
            total_sw_references = sum(len(mw.children) for mw in self.layers[2])
            
            # 冗余度公式: (总引用数 / 唯一节点数 - 1) * 100%
            # 例如: 有 100 个 SW，但 MW 存了 150 个指针，说明有 50 个指针是多重归属的 -> 冗余度 50%
            dag_redundancy = (total_sw_references / sw_count - 1) * 100

        return {
            'raw_count': raw_count,
            'sw_count': sw_count,
            'mw_count': mw_count,
            'lw_count': lw_count,
            
            'compression_sw': comp_sw,
            'compression_mw': comp_mw,
            'compression_lw': comp_lw,
            
            'dag_redundancy': dag_redundancy
        }