"""
搜索引擎
实现 Greedy Search、Exhaustive Search、Adaptive Search 三种策略
"""

from typing import List, Set, Dict, Tuple, Union
from .structure import WorldNode, BioSequence, GenomePointer
from .tools import compute_distance



class SearchStats:
    """搜索性能统计计数器"""
    
    def __init__(self):
        """初始化统计计数器"""
        self.node_access_count = 0
        self.dist_calc_count = 0
        self.leaf_verify_count = 0
        self.layer_breakdown = {
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
        self.index = index_builder
        self.layers = index_builder.layers
    
    def _center_seq(self, node: WorldNode) -> BioSequence:
        """Node 中心转为 BioSequence 以便与 query 计算距离."""
        s = node.get_center_sequence()
        return BioSequence("_center", s)

    def _check_overlap(self, query_seq, node: WorldNode, tolerance: int) -> bool:
        dist = compute_distance(query_seq, self._center_seq(node))
        if dist <= node.radius + tolerance:
            return True, dist
        return False, dist

    def _compute_anchor_dists(self, query_seq: BioSequence, node: WorldNode,
                              stats: SearchStats) -> List[int]:
        """计算 Query 到 node 内所有 routing_anchors 的距离。"""
        dists = []
        for anchor in node.routing_anchors:
            if isinstance(anchor, BioSequence):
                ac = anchor
            else:
                ac = BioSequence("_center", anchor.get_center_sequence())
            d = compute_distance(query_seq, ac)
            stats.dist_calc_count += 1
            dists.append(d)
        return dists

    def _anchor_prunable(self, node: WorldNode, child: WorldNode,
                         tolerance: int, q_anchor_dists: List[int]) -> bool:
        """
        三角不等式安全剪枝：
        |d(Q,A) - d(child.center, A)| 是 d(Q, child.center) 的下界。
        如果任一 anchor 的下界 > child.radius + tolerance，
        则 Query 球与 child 球不可能 overlap，安全剪掉。
        """
        child_id = node._get_child_id(child)
        fingerprint = node.routing_fingerprints.get(child_id)
        if not fingerprint or len(fingerprint) != len(q_anchor_dists):
            return False
        child_r = child.radius
        for i in range(len(q_anchor_dists)):
            if abs(q_anchor_dists[i] - fingerprint[i]) > child_r + tolerance:
                return True
        return False

    # =========================================================================
    # Greedy Search
    # =========================================================================

    def search_greedy(self, query_seq: BioSequence, tolerance: int) -> Tuple[List[BioSequence], SearchStats]:
        """
        贪婪单路径搜索（Candidate Mode）：
        每层选距离最近的一个节点，用锚点剪枝过滤子节点。
        不做叶子精确验证，返回候选集。
        """
        stats = SearchStats()
        current_layer_nodes = list(self.layers[3])
        
        for layer_id in [3, 2, 1]:
            layer_name = {3: 'LW', 2: 'MW', 1: 'SW'}.get(layer_id, 'UNK')
            
            best_node = None
            min_dist = float('inf')
            
            for node in current_layer_nodes:
                d = compute_distance(query_seq, self._center_seq(node))
                stats.dist_calc_count += 1
                stats.node_access_count += 1
                stats.layer_breakdown[layer_name] += 1
                
                if d <= node.radius + tolerance and d < min_dist:
                    min_dist = d
                    best_node = node
            
            if not best_node:
                return [], stats
            
            if layer_id == 1:
                current_layer_nodes = [c for c in best_node.children if isinstance(c, BioSequence)]
                if not current_layer_nodes and best_node.children:
                     current_layer_nodes = list(best_node.children)
                break
            
            q_anchor_dists = self._compute_anchor_dists(query_seq, best_node, stats)
            next_candidates = []
            for child in best_node.children:
                if not isinstance(child, WorldNode):
                    continue
                if not self._anchor_prunable(best_node, child, tolerance, q_anchor_dists):
                    next_candidates.append(child)
            
            current_layer_nodes = next_candidates
        
        return current_layer_nodes, stats

    # =========================================================================
    # Exhaustive Search
    # =========================================================================

    def search_exhaustive(self, query_seq: BioSequence, tolerance: int) -> Tuple[List[BioSequence], SearchStats]:
        """
        穷举搜索：遍历所有满足几何约束的节点，节点级去重。
        中心包含连边下，过滤条件用累积半径补偿以保证完美召回。
        """
        stats = SearchStats()
        unique_results: Dict[str, BioSequence] = {}
        visited_nodes: Set[str] = set()
        
        def traverse(node: WorldNode, current_layer: int):
            if node.node_id in visited_nodes:
                return
            visited_nodes.add(node.node_id)
            
            dist = compute_distance(query_seq, self._center_seq(node))
            stats.dist_calc_count += 1
            stats.node_access_count += 1
            
            layer_name = {3: 'LW', 2: 'MW', 1: 'SW'}.get(current_layer, 'UNK')
            stats.layer_breakdown[layer_name] += 1
            
            if dist > node.radius + tolerance:
                return
            
            for child in node.children:
                if isinstance(child, WorldNode):
                    traverse(child, current_layer - 1)
                elif current_layer == 1:
                    leaf_dist = compute_distance(query_seq, child)
                    stats.dist_calc_count += 1
                    stats.leaf_verify_count += 1
                    if leaf_dist <= tolerance:
                        child_id = getattr(child, 'id', str(id(child)))
                        unique_results[child_id] = child

        for lw_node in self.layers[3]:
            traverse(lw_node, 3)
            
        return list(unique_results.values()), stats

    # =========================================================================
    # Adaptive Search
    # =========================================================================

    def search_adaptive(self, query_seq: BioSequence, tolerance: int) -> Tuple[List[BioSequence], SearchStats]:
        """
        自适应搜索：零 False Negative + 两层加速。

        依赖构建阶段的覆盖性不变量（球重叠连边保证）：
          叶子 s 在节点 N 球内 => 从 N 子树可达 s

        加速 1 — 同层 early-stop：
          扫描同层候选时，如果 Query 球完全在某节点球内
          (d + tolerance <= radius)，则只进该节点，跳过同层其余候选。

        加速 2 — 锚点三角不等式剪枝：
          进入节点后，用 fingerprint 对子节点计算距离下界，
          下界 > child.radius + tolerance 的安全剪掉。
        """
        stats = SearchStats()
        unique_results: Dict[str, BioSequence] = {}
        visited_nodes: Set[str] = set()

        def process_node(node: WorldNode, current_layer: int):
            if current_layer == 1:
                for child in node.children:
                    if isinstance(child, BioSequence):
                        leaf_dist = compute_distance(query_seq, child)
                        stats.dist_calc_count += 1
                        stats.leaf_verify_count += 1
                        if leaf_dist <= tolerance:
                            unique_results[child.id] = child
                return

            world_children = [c for c in node.children if isinstance(c, WorldNode)]
            if not world_children:
                return

            child_layer = current_layer - 1

            num_anchors = len(node.routing_anchors)
            if len(world_children) > num_anchors + 1 and num_anchors > 0:
                q_anchor_dists = self._compute_anchor_dists(query_seq, node, stats)
                surviving = [c for c in world_children
                             if not self._anchor_prunable(
                                 node, c, tolerance, q_anchor_dists)]
            else:
                surviving = world_children

            search_layer(surviving, child_layer)

        def search_layer(candidates: List[WorldNode], layer_id: int):
            layer_name = {3: 'LW', 2: 'MW', 1: 'SW'}.get(layer_id, 'UNK')

            contained_node = None
            overlap_nodes = []

            for node in candidates:
                if node.node_id in visited_nodes:
                    continue

                d = compute_distance(query_seq, self._center_seq(node))
                stats.dist_calc_count += 1
                stats.node_access_count += 1
                stats.layer_breakdown[layer_name] += 1

                if d > node.radius + tolerance:
                    continue

                if d + tolerance <= node.radius:
                    contained_node = node
                    break
                else:
                    overlap_nodes.append(node)

            if contained_node:
                visited_nodes.add(contained_node.node_id)
                process_node(contained_node, layer_id)
            else:
                for node in overlap_nodes:
                    if node.node_id not in visited_nodes:
                        visited_nodes.add(node.node_id)
                        process_node(node, layer_id)

        search_layer(list(self.layers[3]), 3)

        return list(unique_results.values()), stats

    # =========================================================================
    # Brute Force
    # =========================================================================

    def search_brute_force(self, query_seq: BioSequence, tolerance: int,
                           all_sequences: List[BioSequence]) -> Tuple[List[BioSequence], SearchStats]:
        """暴力全量扫描，作为零 False Negative 的 ground truth。"""
        stats = SearchStats()
        results = []
        for seq in all_sequences:
            d = compute_distance(query_seq, seq)
            stats.dist_calc_count += 1
            stats.leaf_verify_count += 1
            if d <= tolerance:
                results.append(seq)
        return results, stats

    # =========================================================================
    # Unified Interface
    # =========================================================================

    def search(self, query_seq: BioSequence, tolerance: int, mode: str = 'adaptive') -> Tuple[List[BioSequence], SearchStats]:
        """统一搜索接口"""
        if mode == 'greedy':
            return self.search_greedy(query_seq, tolerance)
        elif mode == 'exhaustive':
            return self.search_exhaustive(query_seq, tolerance)
        elif mode == 'adaptive':
            return self.search_adaptive(query_seq, tolerance)
        else:
            raise ValueError(f"Unknown search mode: {mode}")
