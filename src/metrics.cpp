#ifndef METRICS_H
#define METRICS_H

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include "levenshtein.hpp"


namespace metrics {

int pos_to_int(const std::string& pos_str) {
    try {
        return std::stoi(pos_str);
    } catch (...) {
        // 遇到无效位置时，返回一个负值作为标记
        return -1; 
    }
}

/**
 * @brief 根据 "区域匹配" 规则计算 TP, FP, FN。
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
                                      int& tp, int& fp, int& fn) {
    
    // 1. 定义位置容忍度
    const int pos_tolerance = max_dist;
    
    // 2. 将所有位置从 string 转换为 int
    std::vector<int> truth_pos;
    for (const auto& s : truth) {
        int pos = pos_to_int(s);
        if (pos >= 0) truth_pos.push_back(pos);
    }
    
    std::vector<int> candidate_pos;
    for (const auto& s : candidates) {
        int pos = pos_to_int(s);
        if (pos >= 0) candidate_pos.push_back(pos);
    }
    
    // 3. 初始化计数器和标记
    tp = 0;
    fn = 0;
    fp = 0;
    
    // 如果没有真实值，所有候选者都是 FP
    if (truth_pos.empty()) {
        fp = candidate_pos.size();
        return;
    }
    // 如果没有候选者，所有真实值都是 FN
    if (candidate_pos.empty()) {
        fn = truth_pos.size();
        return;
    }

    // 标记哪些 candidate 被 "命中" (用于计算 FP)
    std::vector<bool> candidate_hit(candidate_pos.size(), false);

    // 4. 计算 TP 和 FN (遍历所有 truth)
    // 规则：对于每个 truth，只要其容错范围内有 *任意* candidate，就算 TP。
    for (int t_pos : truth_pos) {
        bool is_tp = false;
        
        // 检查是否有任何 candidate 落在 [t_pos - tol, t_pos + tol] 范围内
        for (size_t j = 0; j < candidate_pos.size(); ++j) {
            int c_pos = candidate_pos[j];
            
            if (std::abs(c_pos - t_pos) <= pos_tolerance) {
                is_tp = true;
                // 标记这个 candidate 是 "有用" 的 (被命中)
                candidate_hit[j] = true; 
            }
        }
        
        if (is_tp) {
            tp++; // 这个 truth 找到了匹配
        } else {
            fn++; // 这个 truth 没找到匹配
        }
    }
    
    // 5. 计算 FP (遍历所有 candidate)
    // 规则：如果一个 candidate 没有落在任何 truth 的容错范围内，
    
    for (size_t j = 0; j < candidate_hit.size(); ++j) {
        if (!candidate_hit[j]) {
            fp++; // 这个 candidate 是多余的
        }
    }

    // 注意：在这个模型中，TP + FN = truth_pos.size()，这是由定义决定的。
    // TP + FP = candidate_pos.size() 不一定成立。
}

// TP = truth & candidates
int countTP(const std::vector<std::string>& truth,
                   const std::vector<std::string>& candidates) {
    int tp = 0;
    for (auto &t : truth)
        if (std::find(candidates.begin(), candidates.end(), t) != candidates.end())
            tp++;
    return tp;
}

// FP = candidates - truth
int countFP(const std::vector<std::string>& truth,
                   const std::vector<std::string>& candidates) {
    int fp = 0;
    for (auto &c : candidates)
        if (std::find(truth.begin(), truth.end(), c) == truth.end())
            fp++;
    return fp;
}

// FN = truth - candidates
int countFN(const std::vector<std::string>& truth,
                   const std::vector<std::string>& candidates) {
    int fn = 0;
    for (auto &t : truth)
        if (std::find(candidates.begin(), candidates.end(), t) == candidates.end())
            fn++;
    return fn;
}

std::pair<double, int> evaluateDistances(
    const std::string &query,
    const std::vector<std::string> &candidates,
    const std::string &ref)
{
    if (candidates.empty())
        return {0.0, 0};

    double totalDist = 0.0;
    int maxDist = 0;

    for (const auto &c_str : candidates) {
        int pos = std::stoi(c_str); // 将字符串转为整数索引
        if (pos + query.size() > ref.size()) continue; // 防止越界

        std::string candidate = ref.substr(pos, query.size());

        int d = levenshtein(query, candidate);
        totalDist += d;
        if (d > maxDist) maxDist = d;
    }

    double avgDist = totalDist / static_cast<double>(candidates.size());
    return {avgDist, maxDist};
}


// 打印 report
void report(const std::vector<std::string>& truth,
                   const std::vector<std::string>& candidates) {
    int tp = countTP(truth, candidates);
    int fp = countFP(truth, candidates);
    int fn = countFN(truth, candidates);

    double recall = tp / double(tp + fn);
    double precision = tp / double(tp + fp);

    std::cout << "TP: " << tp << " FP: " << fp << " FN: " << fn
              << " Recall: " << recall << " Precision: " << precision << "\n";
}

} // namespace metrics

#endif
