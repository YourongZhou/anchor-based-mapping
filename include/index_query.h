#ifndef INDEX_QUERY_H
#define INDEX_QUERY_H

#include "index_query.h"
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
                   int max_distance_store);


// 判断 query 是否有任意 anchor 距离小于阈值（利用 anchor_index）
bool hasNearbyAnchor(
    const std::string &query,
    const std::unordered_map<std::string, std::vector<std::pair<int,int>>> &anchor_index,
    int max_dist_query);

// 过滤 queries：只保留在 anchor_index 范围内有相近 anchor 的 query
std::vector<std::string> filter_queries_by_anchor_index(
    const std::vector<std::string> &queries,
    const std::unordered_map<std::string, std::vector<std::pair<int,int>>> &anchor_index,
    int max_dist_query);

#endif