#ifndef METRICS_H
#define METRICS_H

#include <string>
#include <vector>
#include <unordered_set>
#include <iostream>

namespace metrics {

/**
 * @brief 使用位置容错计算 TP, FP, FN。
 *
 * @param truth 真实的匹配位置集合 (字符串)。
 * @param candidates 算法找到的匹配位置集合 (字符串)。
 * @param max_dist 原始的容错范围 (MaxDist)。
 * @param tp 输出：True Positives
 * @param fp 输出：False Positives
 * @param fn 输出：False Negatives
 */
void calculate_position_metrics(const std::vector<std::string>& truth,
                                      const std::vector<std::string>& candidates,
                                      int max_dist,
                                      int& tp, int& fp, int& fn);

/// 计算 True Positive 数量
/// truth: 真正应该被检索到的 reference 序列（ground truth）
/// candidates: 实际 candidate 集合
/// 返回 TP 数量
int countTP(const std::vector<std::string> &truth,
              const std::vector<std::string> &candidates);

/// 计算 False Positive 数量
/// 返回 FP 数量
int countFP(const std::vector<std::string> &truth,
              const std::vector<std::string> &candidates);

/// 计算 False Negative 数量
/// 返回 FN 数量
int countFN(const std::vector<std::string> &truth,
              const std::vector<std::string> &candidates);

std::pair<double, int> evaluateDistances(
    const std::string &query,
    const std::vector<std::string> &candidates,
    const std::string &ref);

/// 打印简单的正确性与效率报告
void report(const std::vector<std::string> &truth,
            const std::vector<std::string> &candidates);

} // namespace metrics

#endif
