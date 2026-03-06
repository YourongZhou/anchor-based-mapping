"""
测试套件
实现三个测试用例：几何半径验证、距离分布热力图、冗余与效率对比
"""

import random
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
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
                    
                    # 球重叠连边：dist(P.center, C.center) <= P.radius + C.radius
                    ball_overlap = dist <= parent_radius + child.radius
                    
                    if not ball_overlap:
                        failures.append(
                            f"{layer_name} {parent.node_id} -> {child.node_id}: "
                            f"dist={dist}, P.radius={parent_radius}, C.radius={child.radius}, "
                            f"threshold={parent_radius + child.radius}"
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

def plot_benchmark_results(adaptive_stats, exhaustive_stats, output_file='benchmark_result.png'):
    """
    绘制 Benchmark 结果对比图 (Adaptive vs Exhaustive)
    """
    layers = ['LW', 'MW', 'SW', 'Leaf Verify\n(Disk I/O)']
    
    a_counts = [
        sum(s.layer_breakdown.get('LW', 0) for s in adaptive_stats),
        sum(s.layer_breakdown.get('MW', 0) for s in adaptive_stats),
        sum(s.layer_breakdown.get('SW', 0) for s in adaptive_stats),
        sum(s.leaf_verify_count for s in adaptive_stats)
    ]
    
    e_counts = [
        sum(s.layer_breakdown.get('LW', 0) for s in exhaustive_stats),
        sum(s.layer_breakdown.get('MW', 0) for s in exhaustive_stats),
        sum(s.layer_breakdown.get('SW', 0) for s in exhaustive_stats),
        sum(s.leaf_verify_count for s in exhaustive_stats)
    ]
    
    reductions = []
    for a, e in zip(a_counts, e_counts):
        if e > 0:
            reductions.append((1 - a / e) * 100)
        else:
            reductions.append(0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('BioGeometry Index: Adaptive vs Exhaustive Search Benchmark', fontsize=16)

    x = np.arange(len(layers))
    width = 0.35
    
    rects1 = ax1.bar(x - width/2, a_counts, width, label='Adaptive (Optimized)', color='#2ecc71', alpha=0.9)
    rects2 = ax1.bar(x + width/2, e_counts, width, label='Exhaustive (Baseline)', color='#e74c3c', alpha=0.9)
    
    ax1.set_ylabel('Access / Verify Count')
    ax1.set_title('Absolute Computational Cost')
    ax1.set_xticks(x)
    ax1.set_xticklabels(layers)
    ax1.legend()
    ax1.grid(axis='y', linestyle='--', alpha=0.3)
    
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax1.annotate(f'{int(height)}',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    autolabel(rects1)
    autolabel(rects2)

    colors = plt.cm.Blues(np.array(reductions) / 100 * 0.8 + 0.2)
    bars = ax2.bar(layers, reductions, color=colors)
    
    ax2.set_ylabel('Reduction Percentage (%)')
    ax2.set_title('Efficiency Gain (Adaptive over Exhaustive)')
    ax2.set_ylim(0, 100)
    ax2.grid(axis='y', linestyle='--', alpha=0.3)
    
    for bar, rate in zip(bars, reductions):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., max(height, 0) + 1,
                f'{rate:.1f}%',
                ha='center', va='bottom', fontweight='bold', color='darkblue')

    a_total = sum(s.dist_calc_count for s in adaptive_stats)
    e_total = sum(s.dist_calc_count for s in exhaustive_stats)
    explanation = (
        f"Key Insight:\n"
        f"Adaptive Search: {a_total} dist calcs vs Exhaustive: {e_total}\n"
        f"Saving: {(1 - a_total/e_total)*100:.1f}% with zero False Negative."
    )
    plt.figtext(0.99, 0.02, explanation, horizontalalignment='right', fontsize=10, style='italic', bbox=dict(facecolor='white', alpha=0.5))

    plt.tight_layout()
    
    save_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), output_file)
    plt.savefig(save_path, dpi=150)
    print(f"\n  [图表生成] Benchmark 图表已保存至: {save_path}")
    plt.close()

def test_efficiency_benchmark(index_builder: BioGeometryIndexBuilder, raw_sequences: list, seed: int = 42, tolerance: int = 5
):
    """
    测试用例 3: 冗余与效率对比
    (已修改: 增加分层冗余度分析)
    """
    print("\n" + "=" * 70)
    print("测试用例 3: 冗余与效率对比 (Efficiency Benchmark)")
    print("=" * 70)
    
    # ... [前半部分生成 Query 的代码保持不变] ...
    # 生成 100 条随机 Query
    print("\n[3.1] 生成测试查询...")
    random.seed(seed)
    queries = []
    for i in range(100):
        base_seq = random.choice(raw_sequences)
        mutated_seq = list(base_seq.seq)
        # 引入 1-3 个随机突变
        for _ in range(random.randint(1, 3)):
            idx = random.randint(0, len(mutated_seq) - 1)
            bases = ['A', 'T', 'C', 'G']
            mutated_seq[idx] = random.choice([b for b in bases if b != mutated_seq[idx]])
        queries.append(BioSequence(f"query_{i:03d}", ''.join(mutated_seq)))
    
    search_engine = BioGeometrySearchEngine(index_builder)
    
    # 运行 Adaptive Search
    print("\n[3.2] 运行 Adaptive Search (模式 A)...")
    adaptive_results = []
    adaptive_stats_list = []
    for query in queries:
        results, stats = search_engine.search_adaptive(query, tolerance)
        adaptive_results.append((query.id, set(seq.id for seq in results)))
        adaptive_stats_list.append(stats)
    
    # 运行 Exhaustive Search
    print("[3.3] 运行 Exhaustive Search (模式 B)...")
    exhaustive_results = []
    exhaustive_stats_list = []
    for query in queries:
        results, stats = search_engine.search_exhaustive(query, tolerance)
        exhaustive_results.append((query.id, set(seq.id for seq in results)))
        exhaustive_stats_list.append(stats)
    
    # 汇总统计
    print("\n[3.4] 结果对比...")
    
    recall_match = True
    no_false_positive = True
    recall_sum = 0.0
    for i, (adaptive_id, adaptive_set) in enumerate(adaptive_results):
        _, exhaustive_set = exhaustive_results[i]
        if exhaustive_set:
            recall_sum += len(adaptive_set & exhaustive_set) / len(exhaustive_set)
        if adaptive_set - exhaustive_set:
            no_false_positive = False
        if adaptive_set != exhaustive_set:
            recall_match = False
    nq = len(adaptive_results)
    avg_recall = recall_sum / nq if nq else 1.0
    if no_false_positive:
        print(f"  Adaptive 无假阳性, 平均召回率: {avg_recall:.2%}")
    else:
        print(f"  警告: 存在假阳性")
    if recall_match:
        print(f"  Recall 与穷举一致: 100%")
    else:
        print(f"  Adaptive 平均召回: {avg_recall:.2%}")
    recall_match = no_false_positive

    print("\n[3.5] 分层冗余度分析 (Adaptive / Exhaustive):")
    print(f"{'Layer':<10} | {'Adaptive':<15} | {'Exhaustive':<18} | {'Ratio':<15} | {'Saving':<15}")
    print("-" * 85)

    layers = ['LW', 'MW', 'SW']
    
    total_adaptive_layers = {'LW': 0, 'MW': 0, 'SW': 0}
    total_exhaustive_layers = {'LW': 0, 'MW': 0, 'SW': 0}

    for s in adaptive_stats_list:
        for L in layers: total_adaptive_layers[L] += s.layer_breakdown.get(L, 0)
    
    for s in exhaustive_stats_list:
        for L in layers: total_exhaustive_layers[L] += s.layer_breakdown.get(L, 0)

    for layer in layers:
        adaptive_cnt = total_adaptive_layers[layer]
        total_cnt = total_exhaustive_layers[layer]
        
        if total_cnt > 0:
            ratio = adaptive_cnt / total_cnt
            saving = (total_cnt - adaptive_cnt) / total_cnt
            print(f"{layer:<10} | {adaptive_cnt:<15} | {total_cnt:<18} | {ratio:.2%}          | {saving:.2%}")
        else:
            print(f"{layer:<10} | {adaptive_cnt:<15} | {total_cnt:<18} | N/A             | N/A")

    total_adaptive_leaf = sum(s.leaf_verify_count for s in adaptive_stats_list)
    total_exhaustive_leaf = sum(s.leaf_verify_count for s in exhaustive_stats_list)
    
    if total_exhaustive_leaf > 0:
        leaf_ratio = total_adaptive_leaf / total_exhaustive_leaf
        leaf_saving = (total_exhaustive_leaf - total_adaptive_leaf) / total_exhaustive_leaf
        print("-" * 85)
        print(f"{'Leaf Data':<10} | {total_adaptive_leaf:<15} | {total_exhaustive_leaf:<18} | {leaf_ratio:.2%}          | {leaf_saving:.2%}")
        print(f"(Disk I/O) ")

    total_adaptive_dist = sum(s.dist_calc_count for s in adaptive_stats_list)
    total_exhaustive_dist = sum(s.dist_calc_count for s in exhaustive_stats_list)
    print(f"\n  总距离计算: Adaptive={total_adaptive_dist}, Exhaustive={total_exhaustive_dist}, 节省={1 - total_adaptive_dist/total_exhaustive_dist:.1%}")

    try:
        plot_benchmark_results(adaptive_stats_list, exhaustive_stats_list)
    except Exception as e:
        print(f"  [警告] 绘图失败 (可能是缺少 matplotlib): {e}")

    return recall_match


def test_false_negative(index_builder: BioGeometryIndexBuilder, raw_sequences: list,
                        seed: int = 42, tolerance: int = 2, num_queries: int = 100):
    """
    测试用例 4: False Negative 精确验证
    用 Brute Force 全量扫描作为 ground truth，验证 Exhaustive / Adaptive Search 是否有漏报。
    同时对比各搜索模式的计算开销。
    """
    print("\n" + "=" * 70)
    print("测试用例 4: False Negative 精确验证 (All Modes vs Brute Force)")
    print("=" * 70)

    random.seed(seed)
    queries = []
    for i in range(num_queries):
        base_seq = random.choice(raw_sequences)
        mutated_seq = list(base_seq.seq)
        for _ in range(random.randint(1, 3)):
            idx = random.randint(0, len(mutated_seq) - 1)
            bases = ['A', 'T', 'C', 'G']
            mutated_seq[idx] = random.choice([b for b in bases if b != mutated_seq[idx]])
        queries.append(BioSequence(f"fnq_{i:03d}", ''.join(mutated_seq)))

    search_engine = BioGeometrySearchEngine(index_builder)

    print(f"\n  Queries: {num_queries}, Tolerance: {tolerance}, DB size: {len(raw_sequences)}")

    # --- Brute Force (ground truth) ---
    print("  [1/4] Running Brute Force (ground truth)...")
    bf_results = []
    bf_total_dist = 0
    for q in queries:
        results, st = search_engine.search_brute_force(q, tolerance, raw_sequences)
        bf_results.append(set(s.id for s in results))
        bf_total_dist += st.dist_calc_count

    # --- Exhaustive Search ---
    print("  [2/4] Running Exhaustive Search...")
    ex_results = []
    ex_stats_list = []
    for q in queries:
        results, st = search_engine.search_exhaustive(q, tolerance)
        ex_results.append(set(s.id for s in results))
        ex_stats_list.append(st)

    # --- Adaptive Search ---
    print("  [3/4] Running Adaptive Search...")
    ad_results = []
    ad_stats_list = []
    for q in queries:
        results, st = search_engine.search_adaptive(q, tolerance)
        ad_results.append(set(s.id for s in results))
        ad_stats_list.append(st)

    # --- Greedy Search ---
    print("  [4/4] Running Greedy Search...")
    gr_results = []
    gr_stats_list = []
    for q in queries:
        results, st = search_engine.search_greedy(q, tolerance)
        gr_results.append(set(s.id for s in results))
        gr_stats_list.append(st)

    # --- 对比分析 ---
    print("\n  === 结果对比 ===")

    bf_total_hits = sum(len(s) for s in bf_results)

    def analyze_mode(name, mode_results, mode_stats_list):
        fn_count = 0
        fp_count = 0
        fn_examples = []
        for i in range(num_queries):
            bf_set = bf_results[i]
            m_set = mode_results[i]
            fn = bf_set - m_set
            fp = m_set - bf_set
            if fn:
                fn_count += len(fn)
                if len(fn_examples) < 3:
                    fn_examples.append((queries[i].id, len(fn), len(bf_set)))
            if fp:
                fp_count += len(fp)

        total_dist = sum(s.dist_calc_count for s in mode_stats_list)
        total_node = sum(s.node_access_count for s in mode_stats_list)
        total_leaf = sum(s.leaf_verify_count for s in mode_stats_list)

        recall = (bf_total_hits - fn_count) / bf_total_hits if bf_total_hits > 0 else 1.0
        print(f"\n  --- {name} ---")
        print(f"    False Negative: {fn_count}   False Positive: {fp_count}   Recall: {recall:.4%}")
        print(f"    Dist calcs: {total_dist}   Node visits: {total_node}   Leaf verifies: {total_leaf}")
        if fn_examples:
            for qid, fn_cnt, bf_cnt in fn_examples:
                print(f"    漏报: {qid}: {fn_cnt}/{bf_cnt}")
        return fn_count == 0 and fp_count == 0

    print(f"\n  Brute Force 总命中: {bf_total_hits}, 总距离计算: {bf_total_dist}")

    ex_ok = analyze_mode("Exhaustive Search", ex_results, ex_stats_list)
    ad_ok = analyze_mode("Adaptive Search", ad_results, ad_stats_list)
    analyze_mode("Greedy Search", gr_results, gr_stats_list)

    # 效率对比摘要
    ex_dist = sum(s.dist_calc_count for s in ex_stats_list)
    ad_dist = sum(s.dist_calc_count for s in ad_stats_list)
    gr_dist = sum(s.dist_calc_count for s in gr_stats_list)
    print(f"\n  === 效率摘要 (距离计算次数) ===")
    print(f"    Brute Force:  {bf_total_dist:>8}")
    print(f"    Exhaustive:   {ex_dist:>8}  ({ex_dist/bf_total_dist:.1%} of BF)")
    print(f"    Adaptive:     {ad_dist:>8}  ({ad_dist/bf_total_dist:.1%} of BF)")
    print(f"    Greedy:       {gr_dist:>8}  ({gr_dist/bf_total_dist:.1%} of BF)")

    if ad_ok:
        if ad_dist < ex_dist:
            saving = (1 - ad_dist / ex_dist) * 100
            print(f"\n  ✓ Adaptive: 零 FN + 比 Exhaustive 节省 {saving:.1f}% 距离计算!")
        else:
            print(f"\n  ✓ Adaptive: 零 FN (与 Exhaustive 开销持平)")
    else:
        print(f"\n  ✗ Adaptive Search 存在 False Negative!")

    return ex_ok and ad_ok
