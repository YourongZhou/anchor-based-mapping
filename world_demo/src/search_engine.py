"""
搜索引擎
实现 Greedy Search 和 Exhaustive Search 两种策略
"""

from typing import List, Set, Dict, Tuple, Union
from .structure import WorldNode, BioSequence, GenomePointer
from .tools import compute_distance


class SearchStats:
    """搜索性能统计计数器"""
    
    def __init__(self):
        """初始化统计计数器"""
        self.node_access_count = 0  # 访问的骨架节点数 (LW/MW/SW)
        self.dist_calc_count = 0    # Edit Distance 计算总次数 (含骨架和叶子)
        self.leaf_verify_count = 0  # 最终验证的叶子序列数
        self.layer_breakdown = {    # 各层访问统计
            'LW': 0,
            'MW': 0,
            'SW': 0
        }
    
    def to_dict(self) -> Dict:
        """转换为字典格式"""
        return {
            'node_access_count': self.node_access_count,
            'dist_calc_count': self.dist_calc_count,
            'leaf_verify_count': self.leaf_verify_count,
            'layer_breakdown': self.layer_breakdown.copy()
        }


class BioGeometrySearchEngine:
    """搜索引擎主类"""
    
    def __init__(self, index_builder):
        """
        初始化搜索引擎
        Args:
            index_builder: BioGeometryIndexBuilder 实例
        """
        self.index = index_builder
        # 引用 layers，方便访问
        self.layers = index_builder.layers
    
    def _check_overlap(self, query_seq, node: WorldNode, tolerance: int) -> bool:
        """
        [Helper] 检查 Query 是否与 Node 的覆盖范围重叠
        逻辑: dist(Q, Center) <= Node.Radius + Tolerance
        """
        # 计算中心点距离
        # node.center_ptr 可能是 BioSequence 或 GenomePointer
        dist = compute_distance(query_seq, node.center_ptr)
        
        # 判定重叠
        if dist <= node.radius + tolerance:
            return True, dist
        return False, dist

    def search_greedy(self, query_seq: BioSequence, tolerance: int) -> Tuple[List[BioSequence], SearchStats]:
        """
        贪婪去重搜索（模式 A）
        利用全局 visited_ids Set 避免重复处理同一节点 (DAG 去重)
        """
        stats = SearchStats()
        # 记录已处理过的 WorldNode ID，防止 DAG 多路径重复访问
        visited_node_ids: Set[str] = set()
        results: List[BioSequence] = []
        
        # 从 Layer 3 (LW) 开始
        # 注意：candidates 存储的是"待检查"的节点
        candidates: List[WorldNode] = list(self.layers[3])
        
        # 标记初始节点为已访问 (防止环，虽然这里应该是 DAG)
        for node in candidates:
            visited_node_ids.add(node.node_id)
            
        # 逐层向下: 3(LW) -> 2(MW) -> 1(SW)
        for layer_id in [3, 2, 1]:
            layer_name = {3: 'LW', 2: 'MW', 1: 'SW'}.get(layer_id, 'UNK')
            next_candidates: List[WorldNode] = []
            
            for node in candidates:
                # 1. 几何剪枝 (Pruning)
                # 计算 Query 到 Node 中心的距离
                dist = compute_distance(query_seq, node.center_ptr)
                stats.dist_calc_count += 1
                stats.node_access_count += 1
                stats.layer_breakdown[layer_name] += 1
                
                # 如果不重叠，剪枝 (Skip)
                if dist > node.radius + tolerance:
                    continue
                
                # 2. 处理子节点 (Drill Down)
                for child in node.children:
                    if isinstance(child, WorldNode):
                        # --- 核心去重逻辑 ---
                        # 只有当这个子节点没被访问过时，才加入下一轮
                        if child.node_id not in visited_node_ids:
                            visited_node_ids.add(child.node_id)
                            next_candidates.append(child)
                            
                    elif layer_id == 1: 
                        # SW 层的子节点是数据 (BioSequence / Pointer)
                        # 这里是 Leaf Verify 阶段
                        # 注意：如果是 Pointer，compute_distance 内部需处理提取逻辑
                        leaf_dist = compute_distance(query_seq, child)
                        stats.dist_calc_count += 1
                        stats.leaf_verify_count += 1
                        
                        if leaf_dist <= tolerance:
                            # 如果是 Pointer，这里可能需要转回 Sequence 对象
                            # 假设 child 就是 BioSequence 或有类似接口
                            results.append(child)
            
            # 进入下一层
            candidates = next_candidates
            if not candidates and layer_id > 1:
                break # 没路了，提前结束
        
        return results, stats
    
    def search_exhaustive(self, query_seq: BioSequence, tolerance: int) -> Tuple[List[BioSequence], SearchStats]:
        """
        穷举路径搜索（模式 B）
        不使用 visited_ids，完全按照树/图的路径遍历
        用于对比测试，展示 DAG 的冗余计算量
        """
        stats = SearchStats()
        # 结果集用 Set 存储 ID 以去重 (因为同一结果可能从多条路找到)
        # 但为了避免最后 O(N) 查表，我们同时存储对象
        unique_results: Dict[str, BioSequence] = {}
        
        def traverse(node: WorldNode, current_layer: int):
            # 1. 几何剪枝
            dist = compute_distance(query_seq, node.center_ptr)
            stats.dist_calc_count += 1
            stats.node_access_count += 1
            
            layer_name = {3: 'LW', 2: 'MW', 1: 'SW'}.get(current_layer, 'UNK')
            stats.layer_breakdown[layer_name] += 1
            
            if dist > node.radius + tolerance:
                return # 剪枝
            
            # 2. 遍历所有子路径 (不做 visited 检查)
            for child in node.children:
                if isinstance(child, WorldNode):
                    traverse(child, current_layer - 1)
                elif current_layer == 1:
                    # Leaf Verify
                    leaf_dist = compute_distance(query_seq, child)
                    stats.dist_calc_count += 1
                    stats.leaf_verify_count += 1
                    
                    if leaf_dist <= tolerance:
                        # 假设 child 有 id 属性
                        child_id = getattr(child, 'id', str(id(child)))
                        unique_results[child_id] = child

        # 从所有顶层节点开始递归
        for lw_node in self.layers[3]:
            traverse(lw_node, 3)
            
        return list(unique_results.values()), stats
    
    def search(self, query_seq: BioSequence, tolerance: int, mode: str = 'greedy') -> Tuple[List[BioSequence], SearchStats]:
        """统一搜索接口"""
        if mode == 'greedy':
            return self.search_greedy(query_seq, tolerance)
        elif mode == 'exhaustive':
            return self.search_exhaustive(query_seq, tolerance)
        else:
            raise ValueError(f"Unknown search mode: {mode}")