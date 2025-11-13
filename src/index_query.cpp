#include "fasta_utils_seqan.hpp"
#include "levenshtein.hpp"

#include <seqan/find.h> 

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <cstdlib>


// Anchor 数据结构加上预计算表
struct AnchorIndex {
    FastaRecord anchor;
    std::unordered_map<int, std::vector<int>> distTable; 
};

// anchor 索引：每个 anchor.seq -> (ref_pos, distance)
std::unordered_map<std::string, std::vector<std::pair<int,int>>>
build_anchor_index(const std::vector<FastaRecord> &anchors,
                   const std::string &ref,
                   int max_distance_store) {
    std::unordered_map<std::string, std::vector<std::pair<int,int>>> anchor_index;
    int ref_len = ref.size();
    if (anchors.empty()) return anchor_index;
    int k = anchors[0].seq.size(); // anchor 长度

    int total_substrings = ref_len - k + 1;

    // 对每个 anchor 并行计算
    #pragma omp parallel for schedule(dynamic)
    for (size_t a = 0; a < anchors.size(); a++) {
        const auto &anchor = anchors[a].seq;
        std::vector<std::pair<int,int>> local_matches;
        local_matches.reserve(total_substrings);

        for (int pos = 0; pos < total_substrings; pos++) {
            std::string ref_sub = ref.substr(pos, k);

            // 带 early exit 的编辑距离
            int dist = levenshtein_with_threshold(anchor, ref_sub, max_distance_store);

            if (dist <= max_distance_store) {
                local_matches.emplace_back(pos, dist);
            }
        }

        // 合并到全局哈希表
        #pragma omp critical
        {
        anchor_index[anchor].insert(anchor_index[anchor].end(),
                                    local_matches.begin(),
                                    local_matches.end());
        }
    }
    return anchor_index;
}


// 判断 query 是否有任意 anchor 距离小于阈值（利用 anchor_index）
bool hasNearbyAnchor(
    const std::string &query,
    const std::unordered_map<std::string, std::vector<std::pair<int,int>>> &anchor_index,
    int max_dist_query)
{
    int k = anchor_index.begin()->first.size();  // anchor 长度（假设非空）
    bool found = false;

    // 遍历 anchor_index 的每个 anchor
    for (const auto &[anchor_seq, matches] : anchor_index) {
        // 如果这个 anchor 在 reference 中的匹配距离都比 query 长度大太多，就可以跳过（可选优化）

        // 直接计算 query 与 anchor_seq 的编辑距离（带 early exit）
        int dist = levenshtein_with_threshold(query, anchor_seq, max_dist_query);
        if (dist < max_dist_query) {
            found = true;
            break;
        }
    }
    return found;
}

// 过滤 queries：只保留在 anchor_index 范围内有相近 anchor 的 query
std::vector<std::string> filter_queries_by_anchor_index(
    const std::vector<std::string> &queries,
    const std::unordered_map<std::string, std::vector<std::pair<int,int>>> &anchor_index,
    int max_dist_query)
{
    std::vector<std::string> filtered;
    filtered.reserve(queries.size());

    for (const auto &q : queries) {
        if (hasNearbyAnchor(q, anchor_index, max_dist_query)) {
            filtered.push_back(q);
        }
    }

    std::cout << "Filtered queries: " << filtered.size()
              << " / " << queries.size()
              << " retained ("
              << (100.0 * filtered.size() / queries.size())
              << "%)\n";
    return filtered;
}