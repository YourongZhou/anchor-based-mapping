#ifndef CANDIDATE_RETRIEVER_H
#define CANDIDATE_RETRIEVER_H

#include <seqan/find.h> 

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <nlohmann/json.hpp>
#include "functions.h"
#include "mtree.h"
#include "mtree_types.h"
#include "data_types.h"


// anchor_index: anchor_seq -> vector of (ref_pos, dist_ref)
CandidateResults retrieveCandidates_anchor(
    const std::string &query,
    const std::unordered_map<std::string, std::vector<std::pair<int,int>>> &anchor_index,
    const seqan::Dna5String &ref_seq,
    int k,              // anchor 长度（仅作安全检查）
    int max_dist_query, // 扩展窗口大小（d +/- max_dist_query）
    bool use_all = true, // 是否选择全部 anchor
    int start_idx = 0,  // 起始比例
    int end_idx   = 20,  // 结束比例
    bool last_random = false,  // 最后一个是否随机选择（当按顺序选择的时候） 
    int anchor_radius = 3,
    size_t* dist_count = nullptr,  // 新增指针参数
    bool require_distance = false
);

#endif