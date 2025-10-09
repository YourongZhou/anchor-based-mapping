// include/levenshtein.hpp
#ifndef LEVENSHTEIN_HPP
#define LEVENSHTEIN_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <climits>

// Classic Levenshtein two-row DP
inline int levenshtein(const std::string &a, const std::string &b) {
    size_t n = a.size(), m = b.size();
    if (n == 0) return (int)m;
    if (m == 0) return (int)n;
    std::vector<int> prev(m + 1), cur(m + 1);
    for (size_t j = 0; j <= m; ++j) prev[j] = (int)j;
    for (size_t i = 1; i <= n; ++i) {
        cur[0] = (int)i;
        for (size_t j = 1; j <= m; ++j) {
            int cost = (a[i-1] == b[j-1]) ? 0 : 1;
            int del = prev[j] + 1;
            int ins = cur[j-1] + 1;
            int rep = prev[j-1] + cost;
            cur[j] = std::min(del, std::min(ins, rep));
        }
        prev.swap(cur);
    }
    return prev[m];
}

// Compute minimal edit distance between pattern p and any substring window of text t.
// If t shorter than p, fallback to levenshtein(p, t).
inline int min_edit_distance_window(const std::string &p, const std::string &t) {
    if (t.size() < p.size()) {
        return levenshtein(p, t);
    }
    int best = INT_MAX;
    for (size_t i = 0; i + p.size() <= t.size(); ++i) {
        int d = levenshtein(p, t.substr(i, p.size()));
        if (d < best) {
            best = d;
            if (best == 0) break;
        }
    }
    return best;
}

// 计算带 early exit 的 Levenshtein
inline int levenshtein_with_threshold(const std::string &s1, const std::string &s2, int max_dist) {
    int n = s1.size();
    int m = s2.size();

    if (abs(n - m) > max_dist) return max_dist + 1; // 长度差直接超限
    std::vector<int> prev(m + 1), curr(m + 1);

    for (int j = 0; j <= m; j++) prev[j] = j;

    for (int i = 1; i <= n; i++) {
        curr[0] = i;
        int row_min = curr[0];
        for (int j = 1; j <= m; j++) {
            int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            curr[j] = std::min({ 
                prev[j] + 1,       // deletion
                curr[j - 1] + 1,   // insertion
                prev[j - 1] + cost // substitution
            });
            row_min = std::min(row_min, curr[j]);
        }
        if (row_min > max_dist) return max_dist + 1; // early exit
        std::swap(prev, curr);
    }
    return prev[m];
}

#endif // LEVENSHTEIN_HPP
