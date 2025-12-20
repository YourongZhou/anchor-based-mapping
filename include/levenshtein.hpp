#ifndef LEVENSHTEIN_HPP
#define LEVENSHTEIN_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <string_view>

// 带阈值的 Levenshtein，支持 string_view 以避免内存拷贝
inline int levenshtein_with_threshold(std::string_view s1, std::string_view s2, int max_dist) {
    int n = (int)s1.size();
    int m = (int)s2.size();

    if (std::abs(n - m) > max_dist) return max_dist + 1;
    if (n < m) return levenshtein_with_threshold(s2, s1, max_dist);

    std::vector<int> prev(m + 1);
    std::vector<int> curr(m + 1);

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
            if (curr[j] < row_min) row_min = curr[j];
        }
        if (row_min > max_dist) return max_dist + 1;
        prev.swap(curr);
    }
    return prev[m];
}

inline int levenshtein(std::string_view a, std::string_view b) {
    return levenshtein_with_threshold(a, b, INT_MAX - 1);
}

// Compute minimal edit distance between pattern p and any substring window of text t.
inline int min_edit_distance_window(std::string_view p, std::string_view t) {
    if (t.size() < p.size()) {
        return levenshtein(p, t);
    }
    int best = INT_MAX;
    for (size_t i = 0; i + p.size() <= t.size(); ++i) {
        int d = levenshtein_with_threshold(p, t.substr(i, p.size()), best - 1);
        if (d < best) {
            best = d;
            if (best == 0) break;
        }
    }
    return best;
}

#endif // LEVENSHTEIN_HPP
