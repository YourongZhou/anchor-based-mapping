"""
测试套件
实现三个测试用例：几何半径验证、距离分布热力图、冗余与效率对比
"""

import random
import sys
import os

# 添加父目录到路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.structure import BioSequence, R_SW, R_MW, R_LW
from src.tools import (
    compute_distance,
    generate_reference_sequence,
    generate_reads_with_mutations
)
from src.index_builder import BioGeometryIndexBuilder
from src.search_engine import BioGeometrySearchEngine


def test_radius_compliance(index_builder: BioGeometryIndexBuilder):
    """
    测试用例 1: 几何半径验证
    验证索引的几何一致性
    """
    print("\n" + "=" * 70)
    print("测试用例 1: 几何半径验证 (Radius Compliance)")
    print("=" * 70)
    
    failures = []
    total_checks = 0
    
    # 1. 验证父节点与子节点的包含关系
    print("\n[1.1] 验证父节点与子节点的包含关系...")
    
    for layer_id in [3, 2, 1]:
        layer_name = {3: 'LW', 2: 'MW', 1: 'SW'}.get(layer_id, 'UNK')
        parent_radius = {3: R_LW, 2: R_MW, 1: R_SW}.get(layer_id, 0)
        
        if layer_id == 3:
            # LW 层：检查与 MW 子节点的关系
            parent_nodes = index_builder.layers[3]
            child_layer = index_builder.layers[2]
        elif layer_id == 2:
            # MW 层：检查与 SW 子节点的关系
            parent_nodes = index_builder.layers[2]
            child_layer = index_builder.layers[1]
        else:
            # SW 层：检查与原始序列的关系
            parent_nodes = index_builder.layers[1]
            child_layer = None
        
        for parent in parent_nodes:
            for child in parent.children:
                total_checks += 1
                
                if isinstance(child, BioSequence):
                    # SW 层：检查与原始序列的距离
                    if layer_id == 1:
                        parent_seq = parent.get_center_sequence()
                        child_seq = child.seq
                        dist = compute_distance(
                            BioSequence("temp_p", parent_seq),
                            BioSequence("temp_c", child_seq)
                        )
                        if dist > R_SW:
                            failures.append(f"SW {parent.node_id}: dist({parent_seq[:10]}..., {child_seq[:10]}...) = {dist} > R_SW={R_SW}")
                else:
                    # 检查 WorldNode 子节点
                    parent_seq = parent.get_center_sequence()
                    child_seq = child.get_center_sequence()
                    dist = compute_distance(
                        BioSequence("temp_p", parent_seq),
                        BioSequence("temp_c", child_seq)
                    )
                    
                    # 严格包含：dist(P.center, C.center) + C.radius <= P.radius
                    strict_containment = dist + child.radius <= parent_radius
                    # 中心包含：dist(P.center, C.center) <= P.radius
                    center_containment = dist <= parent_radius
                    
                    if not strict_containment and not center_containment:
                        failures.append(
                            f"{layer_name} {parent.node_id} -> {child.node_id}: "
                            f"dist={dist}, P.radius={parent_radius}, C.radius={child.radius}"
                        )
    
    print(f"  完成：检查了 {total_checks} 个节点关系")
    if failures:
        print(f"  ✗ 失败：{len(failures)} 个违反几何约束")
        print(f"  前 5 个失败示例：")
        for i, fail in enumerate(failures[:5]):
            print(f"    {i+1}. {fail}")
    else:
        print(f"  ✓ 通过：所有节点关系满足几何约束")
    
    # 2. 随机抽取 SW 节点，验证其叶子序列
    print("\n[1.2] 随机抽取 SW 节点，验证叶子序列...")
    
    sw_nodes = index_builder.layers[1]
    if not sw_nodes:
        print("  警告：没有 SW 节点可验证")
        return len(failures) == 0
    
    sample_size = min(100, len(sw_nodes))
    sampled_sws = random.sample(sw_nodes, sample_size)
    
    leaf_failures = 0
    leaf_checks = 0
    
    for sw_node in sampled_sws:
        sw_center = sw_node.get_center_sequence()
        sw_center_obj = BioSequence("temp_sw", sw_center)
        
        for child in sw_node.children:
            if isinstance(child, BioSequence):
                leaf_checks += 1
                dist = compute_distance(sw_center_obj, child)
                if dist > R_SW:
                    leaf_failures += 1
                    if leaf_failures <= 3:
                        print(f"    失败示例: SW {sw_node.node_id}, dist={dist} > R_SW={R_SW}")
    
    print(f"  完成：检查了 {leaf_checks} 个叶子序列")
    if leaf_failures == 0:
        print(f"  ✓ 通过：所有叶子序列满足 R_SW 约束")
    else:
        print(f"  ✗ 失败：{leaf_failures} 个叶子序列违反 R_SW 约束")
    
    total_failures = len(failures) + leaf_failures
    print(f"\n[总结] 总失败数: {total_failures} / {total_checks + leaf_checks}")
    
    return total_failures == 0


def test_distance_heatmap(index_builder: BioGeometryIndexBuilder):
    """
    测试用例 2: 距离分布热力图
    验证"小世界内聚，大世界分离"的特性
    """
    print("\n" + "=" * 70)
    print("测试用例 2: 距离分布热力图 (Distance Heatmap)")
    print("=" * 70)
    
    sw_nodes = index_builder.layers[1]
    if len(sw_nodes) < 10:
        print("  警告：SW 节点数量不足，跳过测试")
        return
    
    # Intra-class (类内距离)
    print("\n[2.1] Intra-class (类内距离) 分析...")
    sample_sws = random.sample(sw_nodes, min(10, len(sw_nodes)))
    
    intra_distances = []
    for sw_node in sample_sws:
        # 获取该 SW 内部的所有序列
        sequences = [child for child in sw_node.children if isinstance(child, BioSequence)]
        
        if len(sequences) < 2:
            continue
        
        # 计算两两距离
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                dist = compute_distance(sequences[i], sequences[j])
                intra_distances.append(dist)
    
    if intra_distances:
        avg_intra = sum(intra_distances) / len(intra_distances)
        min_intra = min(intra_distances)
        max_intra = max(intra_distances)
        print(f"  类内距离统计:")
        print(f"    平均值: {avg_intra:.2f} (预期 < 10)")
        print(f"    最小值: {min_intra}")
        print(f"    最大值: {max_intra}")
        print(f"    样本数: {len(intra_distances)}")
    else:
        print("  警告：无法计算类内距离（序列数量不足）")
        avg_intra = 0
    
    # Inter-class (类间距离)
    print("\n[2.2] Inter-class (类间距离) 分析...")
    sample_sws_inter = random.sample(sw_nodes, min(50, len(sw_nodes)))
    
    inter_distances = []
    sw_centers = []
    for sw_node in sample_sws_inter:
        center_seq = sw_node.get_center_sequence()
        sw_centers.append(BioSequence(f"sw_{sw_node.node_id}", center_seq))
    
    # 计算中心点两两距离
    for i in range(len(sw_centers)):
        for j in range(i + 1, len(sw_centers)):
            dist = compute_distance(sw_centers[i], sw_centers[j])
            inter_distances.append(dist)
    
    if inter_distances:
        avg_inter = sum(inter_distances) / len(inter_distances)
        min_inter = min(inter_distances)
        max_inter = max(inter_distances)
        print(f"  类间距离统计:")
        print(f"    平均值: {avg_inter:.2f} (预期 > 类内距离)")
        print(f"    最小值: {min_inter}")
        print(f"    最大值: {max_inter}")
        print(f"    样本数: {len(inter_distances)}")
        
        # 对比
        if avg_intra > 0:
            ratio = avg_inter / avg_intra if avg_intra > 0 else 0
            print(f"\n  类间/类内距离比: {ratio:.2f} (预期 > 1.0)")
    else:
        print("  警告：无法计算类间距离")
    
    # 尝试绘制热力图（如果可能）
    print("\n[2.3] 距离矩阵可视化...")
    try:
        import matplotlib
        matplotlib.use('Agg')  # 使用非交互式后端
        import matplotlib.pyplot as plt
        import numpy as np
        
        # 构建距离矩阵（仅使用前 20 个 SW 中心点以节省计算）
        matrix_size = min(20, len(sample_sws_inter))
        distance_matrix = np.zeros((matrix_size, matrix_size))
        
        for i in range(matrix_size):
            for j in range(matrix_size):
                if i == j:
                    distance_matrix[i][j] = 0
                else:
                    dist = compute_distance(sw_centers[i], sw_centers[j])
                    distance_matrix[i][j] = dist
        
        # 绘制热力图
        plt.figure(figsize=(10, 8))
        plt.imshow(distance_matrix, cmap='viridis', aspect='auto')
        plt.colorbar(label='Edit Distance')
        plt.title('Inter-class Distance Heatmap (SW Centers)')
        plt.xlabel('SW Node Index')
        plt.ylabel('SW Node Index')
        
        output_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'distance_heatmap.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  ✓ 热力图已保存到: {output_path}")
    except ImportError:
        print("  [跳过] matplotlib 未安装，跳过可视化")
    except Exception as e:
        print(f"  [警告] 可视化失败: {e}")
    
    return True


def test_efficiency_benchmark(index_builder: BioGeometryIndexBuilder, raw_sequences: list, seed: int = 42):
    """
    测试用例 3: 冗余与效率对比
    证明 Greedy Search 的优越性
    """
    print("\n" + "=" * 70)
    print("测试用例 3: 冗余与效率对比 (Efficiency Benchmark)")
    print("=" * 70)
    
    # 生成 100 条随机 Query
    print("\n[3.1] 生成测试查询...")
    random.seed(seed)
    
    queries = []
    for i in range(100):
        # 从原始数据中随机选择一个序列并引入突变
        base_seq = random.choice(raw_sequences)
        mutated_seq = base_seq.seq
        
        # 引入 1-3 个随机突变
        seq_list = list(mutated_seq)
        num_mutations = random.randint(1, 3)
        bases = ['A', 'T', 'C', 'G']
        
        for _ in range(num_mutations):
            idx = random.randint(0, len(seq_list) - 1)
            seq_list[idx] = random.choice([b for b in bases if b != seq_list[idx]])
        
        query = BioSequence(f"query_{i:03d}", ''.join(seq_list))
        queries.append(query)
    
    print(f"  生成了 {len(queries)} 条查询")
    
    # 初始化搜索引擎
    search_engine = BioGeometrySearchEngine(index_builder)
    tolerance = 5
    
    # 运行 Greedy Search
    print("\n[3.2] 运行 Greedy Search (模式 A)...")
    greedy_results = []
    greedy_stats_list = []
    
    for query in queries:
        results, stats = search_engine.search_greedy(query, tolerance)
        greedy_results.append((query.id, set(seq.id for seq in results)))
        greedy_stats_list.append(stats)
    
    # 汇总 Greedy 统计
    total_greedy_node_access = sum(s.node_access_count for s in greedy_stats_list)
    total_greedy_dist_calc = sum(s.dist_calc_count for s in greedy_stats_list)
    total_greedy_leaf_verify = sum(s.leaf_verify_count for s in greedy_stats_list)
    
    # 运行 Exhaustive Search
    print("[3.3] 运行 Exhaustive Search (模式 B)...")
    exhaustive_results = []
    exhaustive_stats_list = []
    
    for query in queries:
        results, stats = search_engine.search_exhaustive(query, tolerance)
        exhaustive_results.append((query.id, set(seq.id for seq in results)))
        exhaustive_stats_list.append(stats)
    
    # 汇总 Exhaustive 统计
    total_exhaustive_node_access = sum(s.node_access_count for s in exhaustive_stats_list)
    total_exhaustive_dist_calc = sum(s.dist_calc_count for s in exhaustive_stats_list)
    total_exhaustive_leaf_verify = sum(s.leaf_verify_count for s in exhaustive_stats_list)
    
    # 对比结果
    print("\n[3.4] 结果对比...")
    
    # 验证 Recall 一致性
    recall_match = True
    for i, (greedy_id, greedy_set) in enumerate(greedy_results):
        exhaustive_id, exhaustive_set = exhaustive_results[i]
        if greedy_set != exhaustive_set:
            recall_match = False
            print(f"  ✗ 查询 {greedy_id}: 结果不一致")
            print(f"    Greedy: {len(greedy_set)} 条结果")
            print(f"    Exhaustive: {len(exhaustive_set)} 条结果")
            break
    
    if recall_match:
        print(f"  ✓ Recall 一致性: 100% (所有查询结果完全一致)")
    
    # 性能对比
    print("\n[3.5] 性能统计对比:")
    print(f"{'指标':<25} | {'Greedy (A)':<15} | {'Exhaustive (B)':<15} | {'改进率':<10}")
    print("-" * 70)
    
    node_reduction = (1 - total_greedy_node_access / total_exhaustive_node_access) * 100 if total_exhaustive_node_access > 0 else 0
    dist_reduction = (1 - total_greedy_dist_calc / total_exhaustive_dist_calc) * 100 if total_exhaustive_dist_calc > 0 else 0
    leaf_reduction = (1 - total_greedy_leaf_verify / total_exhaustive_leaf_verify) * 100 if total_exhaustive_leaf_verify > 0 else 0
    
    print(f"{'Node Access Count':<25} | {total_greedy_node_access:<15} | {total_exhaustive_node_access:<15} | {node_reduction:.1f}%")
    print(f"{'Dist Calc Count':<25} | {total_greedy_dist_calc:<15} | {total_exhaustive_dist_calc:<15} | {dist_reduction:.1f}%")
    print(f"{'Leaf Verify Count':<25} | {total_greedy_leaf_verify:<15} | {total_exhaustive_leaf_verify:<15} | {leaf_reduction:.1f}%")
    
    # 计算冗余度
    if total_exhaustive_dist_calc > 0:
        redundancy = (total_exhaustive_dist_calc - total_greedy_dist_calc) / total_exhaustive_dist_calc * 100
        print(f"\n  DAG 冗余度: {redundancy:.2f}%")
        print(f"  (Greedy 算法消除了 {redundancy:.2f}% 的冗余计算)")
    
    return recall_match
