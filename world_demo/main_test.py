#!/usr/bin/env python3
"""
BioGeometry Static DAG Index (v3.0) - 主测试入口
整合所有测试用例，提供完整的测试验证流程
"""

import sys
import os
import time
import random

# 添加当前目录到路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.structure import BioSequence
from src.tools import (
    generate_reference_sequence,
    generate_reads_with_mutations
)
from src.index_builder import BioGeometryIndexBuilder
from tests.benchmark import (
    test_radius_compliance,
    test_distance_heatmap,
    test_efficiency_benchmark
)


def main():
    """主函数：运行完整的测试套件"""
    print("=" * 70)
    print("BioGeometry Static DAG Index (v3.0) - 测试验证套件")
    print("=" * 70)
    
    # ===== 1. 测试数据生成 =====
    print("\n[阶段 1] 生成测试数据...")
    random.seed(42)
    
    # 生成一条 100,000 bp 的长 DNA 序列
    print("  生成参考序列 (100,000 bp)...")
    reference_seq = generate_reference_sequence(length=100000, seed=42)
    print(f"  参考序列长度: {len(reference_seq)} bp")
    
    # 生成 Reads：随机截取切片 + 1-5% 突变
    print("  生成 Reads (随机截取 + 5% 突变)...")
    num_reads = 100  # 生成 1000 条 reads
    read_length = 30
    mutation_rate = 0.03  # 3% 突变率
    
    t_start = time.perf_counter()
    raw_sequences = generate_reads_with_mutations(
        reference_seq=reference_seq,
        num_reads=num_reads,
        read_length=read_length,
        mutation_rate=mutation_rate,
        seed=42
    )
    t_gen = time.perf_counter() - t_start
    
    print(f"  生成了 {len(raw_sequences)} 条 Reads，耗时 {t_gen:.4f}s")
    print(f"  示例序列: {raw_sequences[0]}")
    
    # ===== 2. 构建索引 =====
    print("\n[阶段 2] 构建索引...")
    index_builder = BioGeometryIndexBuilder()
    
    t_start = time.perf_counter()
    index_builder.build(raw_sequences)
    t_build = time.perf_counter() - t_start
    
    print(f"  索引构建完成，耗时 {t_build:.4f}s")
    
    # 输出索引统计信息
    stats = index_builder.get_statistics()
    print("\n  索引统计信息:")
    print(f"    Layer 0 (Raw Data): {stats['raw_count']} 条序列")
    print(f"    Layer 1 (SW): {stats['sw_count']} 个节点 (radius={index_builder.radius_config[1]})")
    print(f"    Layer 2 (MW): {stats['mw_count']} 个节点 (radius={index_builder.radius_config[2]})")
    print(f"    Layer 3 (LW): {stats['lw_count']} 个节点 (radius={index_builder.radius_config[3]})")
    print(f"\n  压缩率:")
    print(f"    Layer 0 -> Layer 1: {stats['compression_sw']:.2f}x")
    print(f"    Layer 1 -> Layer 2: {stats['compression_mw']:.2f}x")
    print(f"    Layer 2 -> Layer 3: {stats['compression_lw']:.2f}x")
    if 'dag_redundancy' in stats:
        print(f"\n  DAG 冗余度: {stats['dag_redundancy']:.1f}%")
    
    # ===== 3. 运行测试用例 =====
    print("\n[阶段 3] 运行测试用例...")
    
    test_results = {}
    
    # 测试用例 1: 几何半径验证
    try:
        print("\n>>> 开始测试用例 1...")
        result_1 = test_radius_compliance(index_builder)
        test_results['test_1'] = result_1
        print(f">>> 测试用例 1 {'通过' if result_1 else '失败'}")
    except Exception as e:
        print(f">>> 测试用例 1 执行出错: {e}")
        import traceback
        traceback.print_exc()
        test_results['test_1'] = False
    
    # 测试用例 2: 距离分布热力图
    try:
        print("\n>>> 开始测试用例 2...")
        result_2 = test_distance_heatmap(index_builder)
        test_results['test_2'] = result_2
        print(f">>> 测试用例 2 {'通过' if result_2 else '失败'}")
    except Exception as e:
        print(f">>> 测试用例 2 执行出错: {e}")
        import traceback
        traceback.print_exc()
        test_results['test_2'] = False
    
    # 测试用例 3: 冗余与效率对比
    try:
        print("\n>>> 开始测试用例 3...")
        result_3 = test_efficiency_benchmark(index_builder, raw_sequences, seed=42)
        test_results['test_3'] = result_3
        print(f">>> 测试用例 3 {'通过' if result_3 else '失败'}")
    except Exception as e:
        print(f">>> 测试用例 3 执行出错: {e}")
        import traceback
        traceback.print_exc()
        test_results['test_3'] = False
    
    # ===== 4. 测试总结 =====
    print("\n" + "=" * 70)
    print("测试总结")
    print("=" * 70)
    
    print(f"\n测试用例 1 (几何半径验证): {'✓ 通过' if test_results.get('test_1', False) else '✗ 失败'}")
    print(f"测试用例 2 (距离分布热力图): {'✓ 通过' if test_results.get('test_2', False) else '✗ 失败'}")
    print(f"测试用例 3 (冗余与效率对比): {'✓ 通过' if test_results.get('test_3', False) else '✗ 失败'}")
    
    all_passed = all(test_results.values())
    
    if all_passed:
        print("\n✓ 所有测试用例通过！")
    else:
        print("\n✗ 部分测试用例失败，请检查输出信息")
    
    print("\n" + "=" * 70)
    print("测试完成")
    print("=" * 70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
