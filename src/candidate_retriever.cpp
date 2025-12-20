// src/candidate_retriever.cpp
// Input: anchors.fa, anchor_dist_json_file, refs.fasta
// Output: list of candidate reference IDs (stdout), one per line
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "fasta_utils_seqan.hpp"
#include "levenshtein.hpp"
#include "functions.h"
#include "mtree.h"
#include "mtree_types.h"
#include "data_types.h"

#include <seqan/find.h> 

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <nlohmann/json.hpp>
#include <omp.h>

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
) {
    int qlen = (int)query.size();
    if (qlen < k) { 
        if (dist_count) *dist_count = 0;
        return {};
    }

    // ---------- Step 1: 计算 query 与所有 anchor 的距离 ----------
    std::vector<std::pair<std::string,int>> dist_list;
    dist_list.reserve(anchor_index.size());

    for (const auto &kv : anchor_index) {
        const std::string &anchor_seq = kv.first;
        if ((int)anchor_seq.size() != k) continue;
        int d = min_edit_distance_window(anchor_seq, query);

        // ✅ 只保留距离在 anchor_radius 以内的 anchor
        // if (d <= anchor_radius)
        dist_list.emplace_back(anchor_seq, d);
    }

    if (dist_list.empty()) return {};

    // ---------- Step 2: 排序并选前 20% ----------
    std::sort(dist_list.begin(), dist_list.end(),
              [](auto &a, auto &b){ return a.second < b.second; }); // 小于代表从小到大，大于代表从大到小
    size_t keepN = std::max<size_t>(1, dist_list.size() * 0.2);

    std::unordered_set<std::string> selected;
    if (use_all){
         for (auto &p : dist_list) selected.insert(p.first); //选择全部 anchor
    } else{
        // size_t start_idx = std::min(dist_list.size(), 
        //                             std::max<size_t>(0, dist_list.size() * start_idx));
        // size_t end_idx   = std::min(dist_list.size(), 
        //                             std::max<size_t>(0, dist_list.size() * end_idx));
        std::cout << "\nStart and end indexes: " << start_idx << ", " << end_idx;
        std::cout << "\nlength of dist_list: " << dist_list.size();
        if (!last_random){
            for (size_t i = start_idx; i < end_idx; i++) {
                selected.insert(dist_list[i].first);
            }
        } else {
            // 前面按顺序选
            for (size_t i = start_idx; i < end_idx - 1; i++) {
                selected.insert(dist_list[i].first);
            }

            // 随机选最后一个
            if (end_idx > start_idx) {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<size_t> dis(start_idx, end_idx - 1);
                size_t random_i = dis(gen);
                selected.insert(dist_list[random_i].first);
            }
        }
    }

    if (dist_count) *dist_count = dist_list.size();  // 保存数量

    std::vector<std::unordered_set<int>> candidate_sets;
    candidate_sets.reserve(64);
    std::unordered_set<int> intersection;
    bool first_set = true;

    double anchor_query_dist = 0.0;
    // 对每个 anchor 计算它与 query 的距离 d（最小窗口编辑距离）
    for (const auto &kv : anchor_index) {
        const std::string &anchor_seq = kv.first;
        if (!selected.count(anchor_seq)) continue;  // 跳过不在前20%的 anchor
        const auto &ref_list = kv.second; // vector<pair<ref_pos, dist_ref>>

        // 可选：只处理长度匹配的 anchor（通常 anchor_seq.size() == k）
        if ((int)anchor_seq.size() != k) {
            continue;
        }

        // 1) 计算 anchor 与 query 的距离 d（如果 query 比 anchor 长，用最小窗口距离）
        int d = min_edit_distance_window(anchor_seq, query);
        anchor_query_dist += d;

        // 2) 计算允许的 dist_ref 范围
        int low = d - max_dist_query;
        if (low < 0) low = 0;
        int high = d + max_dist_query;

        // 3) 从 anchor 对应的 ref_list 里筛出 dist_ref 在 [low, high] 的 ref_pos
        std::unordered_set<int> current_set;
        current_set.reserve(16);
        for (const auto &pr : ref_list) {
            int ref_pos = pr.first;
            int dist_ref = pr.second;
            if (dist_ref >= low && dist_ref <= high) {
                current_set.insert(ref_pos);
            }
        }

        // if (!current_set.empty()) candidate_sets.emplace_back(std::move(current_set));

        if (current_set.empty()) continue; // 空集跳过

        if (first_set) {
            // 第一个有效的集合，直接作为交集起点
            intersection = std::move(current_set);
            first_set = false;
        } else {
            // 后续集合：计算当前 intersection 和 current_set 的交集
            std::unordered_set<int> next_intersection;
            next_intersection.reserve(std::min(intersection.size(), current_set.size()));
            
            // 确保总是遍历较小的集合以加速
            const auto& smaller = (intersection.size() < current_set.size()) ? intersection : current_set;
            const auto& larger  = (intersection.size() < current_set.size()) ? current_set : intersection;

            for (int pos : smaller) {
                if (larger.count(pos)) {
                    next_intersection.insert(pos);
                }
            }
            intersection.swap(next_intersection); // 更新交集
            if (intersection.empty()) break;      // 提前退出
        }
    }
    // std::cout << "\nAverage anchor - query distance: " << anchor_query_dist / selected.size() << " " << selected.size();

    if (intersection.empty()) return {};

    // // 4) 交集 —— 先按集合大小排序，先交小集合以加速
    // std::sort(candidate_sets.begin(), candidate_sets.end(),
    //           [](const auto &a, const auto &b){ return a.size() < b.size(); });

    // std::unordered_set<int> intersection = std::move(candidate_sets[0]);
    // for (size_t i = 1; i < candidate_sets.size() && !intersection.empty(); ++i) {
    //     std::unordered_set<int> next;
    //     next.reserve(std::min(intersection.size(), candidate_sets[i].size()));
    //     for (int pos : intersection) {
    //         if (candidate_sets[i].count(pos)) next.insert(pos);
    //     }
    //     intersection.swap(next);
    // }

    std::vector<int> positions;
    positions.reserve(intersection.size());
    for (int pos : intersection) positions.push_back(pos);
    std::sort(positions.begin(), positions.end());
    
    double total_distance = 0.0;
    const size_t ref_len = seqan::length(ref_seq);
    
    for (int pos : positions) {
        #pragma omp critical (GlobalLoggerLock)
        {
            globalLogger.accessCandidate("extend");
        }
        // 检查 pos 是否有效（pos >= 0）以及是否会超出 ref_seq 边界
        if ((size_t)pos + qlen > ref_len || pos < 0) {
            // std::cout << "out of bounds" << std::endl;
            // 如果超出边界，跳过此候选位置
            continue; 
        }

        // 1. 使用 SeqAn infix 提取 Reference 子序列
        seqan::Infix<const seqan::Dna5String>::Type ref_infix = 
            seqan::infix(ref_seq, pos, pos + qlen); 
        
        // 2. 将 SeqAn Infix (Dna5) 转换为 std::string
        std::string ref_segment;
        seqan::reserve(ref_segment, seqan::length(ref_infix));
        for (auto c : ref_infix) {
            ref_segment += (char)c; // Dna5 to char conversion
        }
        int dist = levenshtein(query, ref_segment);
        std::cout << "query: " << query << ", ref: " << ref_segment << "distance: " << dist << std::endl;
        
        // 3. 累加
        total_distance += (double)dist;
    }
    
    double average_distance = total_distance / positions.size();

    // 转成有序 vector 返回（方便后续处理）
    return CandidateResults{positions, average_distance};
}

// crude JSON map parser expecting {"id":num,...}
static std::map<std::string,int> parse_simple_json_map(const std::string &s) {
    std::map<std::string,int> out;
    size_t i = 0, n = s.size();
    while (i < n && s[i] != '{') ++i;
    if (i == n) return out;
    ++i;
    while (i < n) {
        while (i<n && isspace((unsigned char)s[i])) ++i;
        if (i<n && s[i]=='}') break;
        if (i<n && s[i]=='"') {
            ++i;
            std::string id;
            while (i<n && s[i]!='"') id.push_back(s[i++]);
            ++i;
            while (i<n && (s[i]==':' || isspace((unsigned char)s[i]))) ++i;
            std::string num;
            while (i<n && (s[i]=='-' || isdigit((unsigned char)s[i]))) num.push_back(s[i++]);
            if (!num.empty()) {
                int val = std::stoi(num);
                out[id] = val;
            }
            while (i<n && s[i] != ',' && s[i] != '}') ++i;
            if (i<n && s[i]==',') ++i;
        } else ++i;
    }
    return out;
}

int candidate_retriever_main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: candidate_retriever <anchors.fa> <anchor_dist_json_file> <refs.fasta>\n";
        return 1;
    }
    std::string anchors_path = argv[1];
    std::string json_file = argv[2];
    std::string refs = argv[3];

    auto anchors = read_fasta_seqan(anchors_path);
    std::map<std::string,std::string> anchor_seq;
    for (const auto &a : anchors) anchor_seq[a.id] = a.seq;

    // read JSON file
    std::ifstream ifs(json_file);
    if (!ifs) {
        std::cerr << "Cannot open " << json_file << std::endl;
        return 1;
    }
    std::string js; std::ostringstream oss;
    oss << ifs.rdbuf();
    js = oss.str();
    auto amap = parse_simple_json_map(js);
    if (amap.empty()) {
        std::cerr << "No anchors parsed from " << json_file << std::endl;
        return 1;
    }

    // build vector of (anchor_seq, threshold)
    std::vector<std::pair<std::string,int>> queries;
    for (auto &kv : amap) {
        auto it = anchor_seq.find(kv.first);
        if (it == anchor_seq.end()) {
            std::cerr << "Warning: anchor id " << kv.first << " not found in anchors.fa\n";
            continue;
        }
        queries.emplace_back(it->second, kv.second);
    }
    if (queries.empty()) {
        std::cerr << "No valid anchors for query\n";
        return 1;
    }

    auto ref_records = read_fasta_seqan(refs);
    size_t cand_count = 0;
    for (const auto &r : ref_records) {
        bool ok = true;
        for (const auto &q : queries) {
            int d = min_edit_distance_window(q.first, r.seq);
            if (d > q.second) { ok = false; break; }
        }
        if (ok) {
            std::cout << r.id << "\n";
            ++cand_count;
        }
    }
    std::cerr << "Candidates found: " << cand_count << std::endl;
    return 0;
}
