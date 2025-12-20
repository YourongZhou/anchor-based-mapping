#ifndef METRICS_H
#define METRICS_H

#include <string>
#include <vector>
#include <unordered_set>
#include <iostream>

namespace metrics {

/**
 * @brief 使用位置容错计算 TP, FP, FN。
 */
void calculate_position_metrics(const std::vector<int>& truth_pos,
                                      const std::vector<int>& candidate_pos,
                                      int max_dist,
                                      int& tp, int& fp, int& fn);

/**
 * @brief 旧接口兼容版本 (不推荐使用，内部会转换)。
 */
void calculate_position_metrics(const std::vector<std::string>& truth,
                                      const std::vector<std::string>& candidates,
                                      int max_dist,
                                      int& tp, int& fp, int& fn);

std::pair<double, int> evaluateDistances(
    const std::string &query,
    const std::vector<int> &candidates,
    const std::string &ref);

std::pair<double, int> evaluateDistances(
    const std::string &query,
    const std::vector<std::string> &candidates,
    const std::string &ref);

} // namespace metrics

#endif
