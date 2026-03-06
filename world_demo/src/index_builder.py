"""
索引构建器 (NavigaMer v6 - Ball-Overlap Coverage Invariant)

核心不变量 — 覆盖性：
  对任意叶子 s 和任意节点 N（任意层），
  如果 dist(s, N.center) <= N.radius，则从 N 沿子树可达 s。

实现方式（三角不等式推导）：
  - 叶子补插：每个叶子出现在所有 dist(s, SW.center) <= R_SW 的 SW 中
  - 球重叠连边：
      SW→MW: dist(SW.center, MW.center) <= R_MW + R_SW
      MW→LW: dist(MW.center, LW.center) <= R_LW + R_MW
    这保证：若叶子 s 在 MW 球内，则包含 s 的某个 SW 一定被连到该 MW。
    证明：dist(SW.center, MW.center) <= dist(SW.center, s) + dist(s, MW.center)
                                      <= R_SW + R_MW
    同理对 MW→LW 层。

搜索时过滤条件：
  由于球重叠连边已保证覆盖性，搜索时只需 dist > N.radius + tolerance 即可安全剪枝，
  无需额外的累积半径补偿。
"""

import random
from typing import List, Dict, Union, Set
from .structure import WorldNode, BioSequence, GenomePointer, R_SW, R_MW, R_LW
from .tools import compute_distance, farthest_point_sampling


def _to_seq_obj(obj):
    """Get an object with .seq for compute_distance. BioSequence -> self; WorldNode -> center as BioSequence."""
    if isinstance(obj, BioSequence):
        return obj
    if isinstance(obj, WorldNode):
        return BioSequence("_center", obj.get_center_sequence())
    if isinstance(obj, GenomePointer):
        return BioSequence("_center", obj.get_sequence())
    return None


class BioGeometryIndexBuilder:
    """索引构建器主类 (NavigaMer 穷举插入 + 路由锚点)"""
    
    ANCHOR_COUNT = 3  # 每节点内部锚点数量
    
    def __init__(self):
        """初始化索引构建器"""
        self.layers: Dict[int, List[WorldNode]] = {
            1: [],
            2: [],
            3: []
        }
        self.radius_config = {1: R_SW, 2: R_MW, 3: R_LW}
        self.stats = {
            'added_sequences': 0,
            'created_nodes': {1: 0, 2: 0, 3: 0}
        }
    
    def build(self, raw_sequences: List[BioSequence]):
        """
        三阶段构建：
        Phase 1 - 骨架构建（增量递归插入，创建节点）
        Phase 2 - 全局 DAG 连通性优化（穷举重连，保证几何包含即连接）
        Phase 3 - 路由锚点精炼（FPS 选锚点 + 重算指纹）
        """
        print(f"[Build] Starting NavigaMer Build for {len(raw_sequences)} sequences...")
        
        # --- Phase 1: Skeleton Build ---
        shuffled_seqs = list(raw_sequences)
        random.shuffle(shuffled_seqs)
        for i, seq in enumerate(shuffled_seqs):
            self.stats['added_sequences'] += 1
            self._recursive_insert(seq, self.layers[3])
            if (i + 1) % 100 == 0:
                print(f"  Phase 1: Processed {i + 1} sequences...", end='\r')
        print(f"\n  Phase 1 done. SW={len(self.layers[1])}, MW={len(self.layers[2])}, LW={len(self.layers[3])}")
        
        # --- Phase 2: Global DAG Connectivity Optimization ---
        print("  Phase 2: Optimizing DAG connectivity (Exhaustive Relinking)...")
        self.optimize_dag_connectivity()
        
        # --- Phase 3: Anchor Refinement ---
        print("  Phase 3: Refining routing anchors...")
        self.refine_index()
        
        print(f"[Build] Completed. Added {self.stats['added_sequences']} sequences.")
        self._print_summary()

    def _recursive_insert(self, item: BioSequence, candidates: List[WorldNode]):
        """
        穷举递归：找到所有包含 item 的父节点，在目标层向所有有效父节点插入；否则创建新节点并向上冒泡。
        """
        item_center = item  # BioSequence 自身即中心
        item_radius = 0
        valid_parents = []
        for node in candidates:
            q = _to_seq_obj(node.center_ptr)
            if q is None:
                continue
            d = compute_distance(item_center, q)
            if d <= node.radius:
                valid_parents.append(node)
        
        if not valid_parents:
            self._handle_new_node_creation(item)
            return
        
        current_layer = valid_parents[0].layer
        target_layer = 1  # 序列插入目标为 SW 层
        
        if current_layer == target_layer:
            for parent in valid_parents:
                self._add_child_to_node(parent, item)
        else:
            next_candidates = []
            seen = set()
            for parent in valid_parents:
                for c in parent.children:
                    if isinstance(c, WorldNode):
                        if c.node_id not in seen:
                            seen.add(c.node_id)
                            next_candidates.append(c)
            self._recursive_insert(item, next_candidates)

    def _handle_new_node_creation(self, item: BioSequence):
        """无父节点时：创建新 SW，并递归向上插入."""
        new_sw = WorldNode(center_ptr=item, radius=R_SW, layer_level=1)
        self._add_child_to_node(new_sw, item)
        new_sw.data_count = 1
        self.layers[1].append(new_sw)
        self.stats['created_nodes'][1] += 1
        self._insert_node_upwards(new_sw)

    def _add_child_to_node(self, parent: WorldNode, child: Union[WorldNode, BioSequence]):
        """带锚点与指纹的子节点添加（增量前 K 锚点策略）."""
        if len(parent.routing_anchors) < self.ANCHOR_COUNT:
            parent.routing_anchors.append(child)
        
        child_center = child if isinstance(child, BioSequence) else _to_seq_obj(child)
        if child_center is None:
            child_center = _to_seq_obj(getattr(child, 'center_ptr', child))
        fingerprint = []
        for anchor in parent.routing_anchors:
            anchor_center = anchor if isinstance(anchor, BioSequence) else _to_seq_obj(anchor)
            if anchor_center is None:
                anchor_center = _to_seq_obj(getattr(anchor, 'center_ptr', anchor))
            if anchor is child:
                fingerprint.append(0)
            else:
                fingerprint.append(compute_distance(child_center, anchor_center))
        
        parent.add_child_with_fingerprint(child, fingerprint)
        if isinstance(child, BioSequence):
            parent.data_count = len([c for c in parent.children if isinstance(c, BioSequence)])

    def _insert_node_upwards(self, child_node: WorldNode):
        """递归向上插入新节点；找到的父节点用 _add_child_to_node 连接."""
        current_layer = child_node.layer
        if current_layer == 3:
            self.layers[3].append(child_node)
            self.stats['created_nodes'][3] += 1
            return
        
        target_parent_layer = current_layer + 1
        parent_radius = self.radius_config[target_parent_layer]
        valid_parents = self._search_candidates(
            query_center=child_node.center_ptr,
            query_radius=child_node.radius,
            target_layer=target_parent_layer
        )
        
        if valid_parents:
            for parent in valid_parents:
                self._add_child_to_node(parent, child_node)
            return
        
        new_parent = WorldNode(
            center_ptr=child_node.center_ptr,
            radius=parent_radius,
            layer_level=target_parent_layer
        )
        self._add_child_to_node(new_parent, child_node)
        self.layers[target_parent_layer].append(new_parent)
        self.stats['created_nodes'][target_parent_layer] += 1
        self._insert_node_upwards(new_parent)

    def _search_candidates(self, query_center, query_radius: int, target_layer: int) -> List[WorldNode]:
        """自顶向下搜索严格包含 query 的目标层节点（用于 _insert_node_upwards）."""
        current_candidates = self.layers[3]
        current_layer_level = 3
        q_obj = _to_seq_obj(query_center)
        if q_obj is None:
            return []
        
        while current_layer_level > target_layer:
            next_candidates = []
            for node in current_candidates:
                d = compute_distance(_to_seq_obj(node.center_ptr), q_obj)
                if d <= node.radius + query_radius:
                    next_candidates.extend([c for c in node.children if isinstance(c, WorldNode)])
            unique = {n.node_id: n for n in next_candidates}
            current_candidates = list(unique.values())
            current_layer_level -= 1
            if not current_candidates:
                return []
        
        final = []
        for node in current_candidates:
            d = compute_distance(_to_seq_obj(node.center_ptr), q_obj)
            if d + query_radius <= node.radius:
                final.append(node)
        return final

    def optimize_dag_connectivity(self):
        """
        Phase 2: 穷举重连 DAG，建立覆盖性不变量。

        Step A: 叶子补插 — 对每个叶子 s 穷举所有 SW，
                dist(s, SW.center) <= R_SW 则插入。
        Step B: SW -> MW 球重叠连边: dist(SW.center, MW.center) <= R_MW + R_SW
        Step C: MW -> LW 球重叠连边: dist(MW.center, LW.center) <= R_LW + R_MW
        Step D: 清理空节点
        """
        sw_nodes = self.layers[1]
        mw_nodes = self.layers[2]
        lw_nodes = self.layers[3]

        # --- A. 叶子补插：保证每个叶子出现在所有应包含它的 SW 中 ---
        all_leaves: Dict[str, BioSequence] = {}
        for sw in sw_nodes:
            for child in sw.children:
                if isinstance(child, BioSequence):
                    all_leaves[child.id] = child

        print(f"    Leaf backfill: {len(all_leaves)} leaves x {len(sw_nodes)} SW...")
        backfill_count = 0
        for sw in sw_nodes:
            existing_ids = set(
                c.id for c in sw.children if isinstance(c, BioSequence)
            )
            sw_center = _to_seq_obj(sw.center_ptr)
            if sw_center is None:
                continue
            for leaf_id, leaf in all_leaves.items():
                if leaf_id in existing_ids:
                    continue
                d = compute_distance(leaf, sw_center)
                if d <= sw.radius:
                    sw.children.append(leaf)
                    existing_ids.add(leaf_id)
                    backfill_count += 1
            sw.data_count = len([c for c in sw.children if isinstance(c, BioSequence)])
        print(f"    Backfilled {backfill_count} leaf-SW links.")

        # --- B. 重连 SW -> MW (球重叠: dist(SW.center, MW.center) <= R_MW + R_SW) ---
        print(f"    Relinking {len(sw_nodes)} SW -> {len(mw_nodes)} MW (ball overlap)...")
        for mw in mw_nodes:
            mw.children = [c for c in mw.children if isinstance(c, BioSequence)]
        for sw in sw_nodes:
            sw_center = _to_seq_obj(sw.center_ptr)
            if sw_center is None:
                continue
            for mw in mw_nodes:
                mw_center = _to_seq_obj(mw.center_ptr)
                if mw_center is None:
                    continue
                d = compute_distance(sw_center, mw_center)
                if d <= mw.radius + sw.radius:
                    mw.children.append(sw)

        # --- C. 重连 MW -> LW (球重叠: dist(MW.center, LW.center) <= R_LW + R_MW) ---
        print(f"    Relinking {len(mw_nodes)} MW -> {len(lw_nodes)} LW (ball overlap)...")
        for lw in lw_nodes:
            lw.children = []
        for mw in mw_nodes:
            if not mw.children:
                continue
            mw_center = _to_seq_obj(mw.center_ptr)
            if mw_center is None:
                continue
            for lw in lw_nodes:
                lw_center = _to_seq_obj(lw.center_ptr)
                if lw_center is None:
                    continue
                d = compute_distance(mw_center, lw_center)
                if d <= lw.radius + mw.radius:
                    lw.children.append(mw)

        # --- D. 清理空节点 ---
        self.layers[2] = [n for n in mw_nodes if n.children]
        self.layers[3] = [n for n in lw_nodes if n.children]
        print(f"    Done. MW={len(self.layers[2])}, LW={len(self.layers[3])}")

    def refine_index(self):
        """后置精炼：用 FPS 重选锚点并重算所有子节点指纹."""
        all_nodes = []
        for layer in self.layers.values():
            all_nodes.extend(layer)
        
        for node in all_nodes:
            if len(node.children) <= self.ANCHOR_COUNT:
                node.routing_anchors = list(node.children)
            else:
                node.routing_anchors = farthest_point_sampling(
                    node.children, self.ANCHOR_COUNT, compute_distance
                )
            node.routing_fingerprints = {}
            for child in node.children:
                child_center = child if isinstance(child, BioSequence) else _to_seq_obj(child)
                if child_center is None:
                    child_center = _to_seq_obj(getattr(child, 'center_ptr', child))
                fp = []
                for anchor in node.routing_anchors:
                    anchor_center = anchor if isinstance(anchor, BioSequence) else _to_seq_obj(anchor)
                    if anchor_center is None:
                        anchor_center = _to_seq_obj(getattr(anchor, 'center_ptr', anchor))
                    if anchor is child:
                        fp.append(0)
                    else:
                        fp.append(compute_distance(child_center, anchor_center))
                node.routing_fingerprints[node._get_child_id(child)] = fp

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