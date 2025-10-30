#include "tools.h"
#include "levenshtein.hpp"
#include "rng.h"

#include <seqan/find.h> 

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <nlohmann/json.hpp>

// 从 reference 随机生成 num_queries 条长度为 query_len 的子串（模拟 reads）
std::vector<std::string> simulate_queries(const std::string &ref,
                                                 int num_queries,
                                                 size_t query_len) {
    if (ref.size() < query_len) {
        throw std::runtime_error("Reference shorter than query length!");
    }

    std::vector<std::string> queries;
    queries.reserve(num_queries);

    // static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, ref.size() - query_len);

    for (int i = 0; i < num_queries; ++i) {
        size_t start = dist(rng);
        queries.push_back(ref.substr(start, query_len));
    }
    return queries;
}

// 用 SeqAn2 找到 query 在 ref 中的所有 exact-match 位置
std::vector<int> find_all_occurrences(const std::string &query,
                                             const std::string &ref) {
    using namespace seqan;
    typedef String<char> TText;
    typedef String<char> TPattern;

    TText text = ref;
    TPattern pattern = query;

    Finder<TText> finder(text);
    Pattern<TPattern, Horspool> patt(pattern);

    std::vector<int> positions;
    while (find(finder, patt)) {
        positions.push_back(position(finder));
    }
    return positions;
}

// 找到所有 query 在某个距离以内的 neighbor
std::vector<int> find_all_occurrences_approx(
    const std::string &query,
    const std::string &ref,
    int max_dist
) {
    std::vector<int> positions;
    int qlen = query.size();
    int rlen = ref.size();

    if (qlen > rlen) return positions;

    for (int i = 0; i <= rlen - qlen; i++) {
        std::string window = ref.substr(i, qlen);
        int d = levenshtein(query, window);
        if (d <= max_dist) {
            positions.push_back(i);
        }
    }

    return positions;
}

// 小工具：把 int positions 转成 string id （用于复用现有 string-based metrics）
std::vector<std::string> posVecToStrVec(const std::vector<int> &pos) {
    std::vector<std::string> out;
    out.reserve(pos.size());
    for (int p : pos) out.push_back(std::to_string(p));
    return out;
}
