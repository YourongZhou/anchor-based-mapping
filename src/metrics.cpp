#include "metrics.h"
#include "levenshtein.hpp"
#include <algorithm>
#include <cmath>

namespace metrics {

void calculate_position_metrics(const std::vector<int>& truth_pos_unsorted,
                                      const std::vector<int>& candidate_pos_unsorted,
                                      int max_dist,
                                      int& tp, int& fp, int& fn) {
    
    const int pos_tolerance = max_dist;
    
    tp = 0;
    fn = 0;
    fp = 0;
    
    if (truth_pos_unsorted.empty()) {
        fp = (int)candidate_pos_unsorted.size();
        return;
    }
    if (candidate_pos_unsorted.empty()) {
        fn = (int)truth_pos_unsorted.size();
        return;
    }

    std::vector<int> truth_pos = truth_pos_unsorted;
    std::sort(truth_pos.begin(), truth_pos.end());
    
    std::vector<int> candidate_pos = candidate_pos_unsorted;
    std::sort(candidate_pos.begin(), candidate_pos.end());

    std::vector<bool> candidate_hit(candidate_pos.size(), false);

    // 计算 TP 和 FN (对每个 truth，找最近的 candidate)
    for (int t_pos : truth_pos) {
        bool found = false;
        // 使用 binary search 找到第一个可能匹配的 candidate
        auto it = std::lower_bound(candidate_pos.begin(), candidate_pos.end(), t_pos - pos_tolerance);
        
        while (it != candidate_pos.end() && *it <= t_pos + pos_tolerance) {
            found = true;
            candidate_hit[std::distance(candidate_pos.begin(), it)] = true;
            ++it;
        }
        
        if (found) tp++;
        else fn++;
        }
    
    for (bool hit : candidate_hit) {
        if (!hit) fp++;
        }
    }

// 兼容版本
void calculate_position_metrics(const std::vector<std::string>& truth,
                                      const std::vector<std::string>& candidates,
                                      int max_dist,
                                      int& tp, int& fp, int& fn) {
    std::vector<int> t_pos, c_pos;
    for (const auto& s : truth) { try { t_pos.push_back(std::stoi(s)); } catch(...) {} }
    for (const auto& s : candidates) { try { c_pos.push_back(std::stoi(s)); } catch(...) {} }
    calculate_position_metrics(t_pos, c_pos, max_dist, tp, fp, fn);
}

std::pair<double, int> evaluateDistances(
    const std::string &query,
    const std::vector<int> &candidates,
    const std::string &ref)
{
    if (candidates.empty())
        return {0.0, 0};

    double totalDist = 0.0;
    int maxD = 0;
    int count = 0;

    for (int pos : candidates) {
        if (pos < 0 || (size_t)pos + query.size() > ref.size()) continue;

        std::string candidate = ref.substr(pos, query.size());
        int d = levenshtein(query, candidate);
        totalDist += d;
        if (d > maxD) maxD = d;
        count++;
    }

    double avgDist = (count > 0) ? (totalDist / count) : 0.0;
    return {avgDist, maxD};
}

std::pair<double, int> evaluateDistances(
    const std::string &query,
    const std::vector<std::string> &candidates,
    const std::string &ref)
{
    std::vector<int> c_pos;
    for (const auto& s : candidates) { try { c_pos.push_back(std::stoi(s)); } catch(...) {} }
    return evaluateDistances(query, c_pos, ref);
}

} // namespace metrics
