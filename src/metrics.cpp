#ifndef METRICS_H
#define METRICS_H

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include "levenshtein.hpp"


namespace metrics {

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
    const std::vector<std::string> &candidates)
{
    if (candidates.empty())
        return {0.0, 0};

    double totalDist = 0.0;
    int maxDist = 0;

    for (const auto &c : candidates) {
        int d = levenshtein(query, c);
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
