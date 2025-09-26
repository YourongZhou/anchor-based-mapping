#ifndef METRICS_H
#define METRICS_H

#include <string>
#include <vector>
#include <unordered_set>
#include <iostream>

namespace metrics {

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

/// 打印简单的正确性与效率报告
void report(const std::vector<std::string> &truth,
            const std::vector<std::string> &candidates);

} // namespace metrics

#endif
