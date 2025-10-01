// main.cpp

#include "fasta_utils_seqan.hpp"
#include "fastq_utils_seqan.hpp"
#include "levenshtein.hpp"
#include "metrics.h"
#include "anchor_gen.h"

// #define SEQAN_NO_INCLUDE_OMP
#include <seqan/find.h> 

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <random>

// -----------------------------
// helpers
// -----------------------------

// 从 reference 随机生成 num_queries 条长度为 query_len 的子串（模拟 reads）
inline std::vector<std::string> simulate_queries(const std::string &ref,
                                                 int num_queries,
                                                 size_t query_len) {
    if (ref.size() < query_len) {
        throw std::runtime_error("Reference shorter than query length!");
    }

    std::vector<std::string> queries;
    queries.reserve(num_queries);

    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, ref.size() - query_len);

    for (int i = 0; i < num_queries; ++i) {
        size_t start = dist(rng);
        queries.push_back(ref.substr(start, query_len));
    }
    return queries;
}

// 用 SeqAn2 找到 query 在 ref 中的所有 exact-match 位置
inline std::vector<int> find_all_occurrences(const std::string &query,
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

// 思路：计算 query 与每个 anchor 的编辑距离，若 <= maxDist，则把该 anchor 在 ref 中出现的所有位置加入候选集合。
inline std::vector<int> retrieveCandidates(const std::string &query,
                                           const std::vector<FastaRecord> &anchors,
                                           const std::string &ref,
                                           int maxDist = 5) {
    std::unordered_set<int> cand_set;
    for (const auto &a : anchors) {
        int d = levenshtein(a.seq, query);
        if (d <= maxDist) {
            // 找到这个 anchor 在 ref 中出现的所有位置
            auto pos_list = find_all_occurrences(a.seq, ref);
            for (int p : pos_list) cand_set.insert(p);
        }
    }
    // to vector + sort
    std::vector<int> cands;
    cands.reserve(cand_set.size());
    for (int p : cand_set) cands.push_back(p);
    std::sort(cands.begin(), cands.end());
    return cands;
}

// Anchor 数据结构加上预计算表
struct AnchorIndex {
    FastaRecord anchor;
    std::unordered_map<int, std::vector<int>> distTable; 
};

// #include <vector>
// #include <string>
// #include <unordered_map>
// #include <limits>
// #include <omp.h>
// #include <seqan/align.h>

// 计算带 early exit 的 Levenshtein
int levenshtein_with_threshold(const std::string &s1, const std::string &s2, int max_dist) {
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
    // #pragma omp parallel for schedule(dynamic)
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
        // #pragma omp critical
        // {
        anchor_index[anchor].insert(anchor_index[anchor].end(),
                                    local_matches.begin(),
                                    local_matches.end());
        // }
    }
    return anchor_index;
}



// anchor_index: anchor_seq -> vector of (ref_pos, dist_ref)
std::vector<int> retrieveCandidates_anchor(
    const std::string &query,
    const std::unordered_map<std::string, std::vector<std::pair<int,int>>> &anchor_index,
    int k,              // anchor 长度（仅作安全检查）
    int max_dist_query  // 扩展窗口大小（d +/- max_dist_query）
) {
    int qlen = (int)query.size();
    if (qlen < k) return {};

    std::vector<std::unordered_set<int>> candidate_sets;
    candidate_sets.reserve(64);

    // 对每个 anchor 计算它与 query 的距离 d（最小窗口编辑距离）
    for (const auto &kv : anchor_index) {
        const std::string &anchor_seq = kv.first;
        const auto &ref_list = kv.second; // vector<pair<ref_pos, dist_ref>>

        // 可选：只处理长度匹配的 anchor（通常 anchor_seq.size() == k）
        if ((int)anchor_seq.size() != k) {
            continue;
        }

        // 1) 计算 anchor 与 query 的距离 d（如果 query 比 anchor 长，用最小窗口距离）
        int d = min_edit_distance_window(anchor_seq, query);

        // 2) 计算允许的 dist_ref 范围
        int low = d - max_dist_query;
        if (low < 0) low = 0;
        int high = d + max_dist_query;

        // 3) 从 anchor 对应的 ref_list 里筛出 dist_ref 在 [low, high] 的 ref_pos
        std::unordered_set<int> s;
        s.reserve(16);
        for (const auto &pr : ref_list) {
            int ref_pos = pr.first;
            int dist_ref = pr.second;
            if (dist_ref >= low && dist_ref <= high) {
                s.insert(ref_pos);
            }
        }

        if (!s.empty()) candidate_sets.emplace_back(std::move(s));
    }

    if (candidate_sets.empty()) return {};

    // 4) 交集 —— 先按集合大小排序，先交小集合以加速
    std::sort(candidate_sets.begin(), candidate_sets.end(),
              [](const auto &a, const auto &b){ return a.size() < b.size(); });

    std::unordered_set<int> intersection = std::move(candidate_sets[0]);
    for (size_t i = 1; i < candidate_sets.size() && !intersection.empty(); ++i) {
        std::unordered_set<int> next;
        next.reserve(std::min(intersection.size(), candidate_sets[i].size()));
        for (int pos : intersection) {
            if (candidate_sets[i].count(pos)) next.insert(pos);
        }
        intersection.swap(next);
    }

    // 转成有序 vector 返回（方便后续处理）
    std::vector<int> result;
    result.reserve(intersection.size());
    for (int pos : intersection) result.push_back(pos);
    std::sort(result.begin(), result.end());
    return result;
}


// 小工具：把 int positions 转成 string id （用于复用现有 string-based metrics）
inline std::vector<std::string> posVecToStrVec(const std::vector<int> &pos) {
    std::vector<std::string> out;
    out.reserve(pos.size());
    for (int p : pos) out.push_back(std::to_string(p));
    return out;
}

int main() {
    // ----- 参数 -----
    const std::string fasta_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fasta/ecoli.fa";
    const size_t truncate_ref_len = 10000;   // 只取前 10000 个碱基
    const size_t anchor_len = 20;            // anchor 长度
    const int num_anchors = 100;           // anchor 个数
    const int num_queries = 100;             // query 数量
    const size_t query_len = 20;             // query 长度
    const int maxDist = 2;                   // 容差距离

    // ----- 读取 FASTA 并截断 -----
    std::vector<FastaRecord> records;
    try {
        records = read_fasta_seqan(fasta_path);
    } catch (std::runtime_error &e) {
        std::cerr << "Error reading FASTA: " << e.what() << "\n";
        return 1;
    }
    if (records.empty()) {
        std::cerr << "No sequences in FASTA\n";
        return 1;
    }
    std::string ref = records[0].seq.substr(0, std::min(records[0].seq.size(), truncate_ref_len));
    std::cout << "Using reference (truncated) length = " << ref.size() << "\n";

    // ----- 生成 anchors -----
    auto anchors = generate_anchors_from_seq(ref, "ref1", anchor_len, num_anchors);
    std::cout << "Generated " << anchors.size() << " anchors\n";

    // ----- 生成 queries -----
    auto queries = simulate_queries(ref, num_queries, query_len);
    std::cout << "Simulated " << queries.size() << " queries\n";

    // ----- ground truth -----
    std::cout << "Generating ground truths:\n";
    std::vector<std::vector<int>> truth_positions;
    truth_positions.reserve(queries.size());
    for (const auto &q : queries) {
        truth_positions.push_back(find_all_occurrences(q, ref));
    }

    // ----- 构建 anchor index -----
    std::cout << "Building anchor index:\n";
    auto anchor_index = build_anchor_index(anchors, ref, anchor_len);

    // ----- 查询 + 评估 -----
    int totalTP = 0, totalFP = 0, totalFN = 0;
    for (size_t i = 0; i < queries.size(); ++i) {
        const auto &q = queries[i];
        const auto &truth_pos = truth_positions[i];

        // === 用 anchor-based 检索 ===
        auto candidates_pos = retrieveCandidates_anchor(q, anchor_index, anchors[0].seq.size(), maxDist);

        // === 打印 ===
        std::cout << "\n=== Query [" << i << "] ===\n";
        std::cout << "Query seq: " << q << "\n";

        std::cout << "Truth positions (" << truth_pos.size() << "): ";
        for (auto p : truth_pos) std::cout << p << " ";
        std::cout << "\n";
        for (auto p : truth_pos) {
            if (p + q.size() <= ref.size()) {
                std::cout << "  Truth seq @ " << p << ": "
                          << ref.substr(p, q.size()) << "\n";
            }
        }

        std::cout << "Candidate positions (" << candidates_pos.size() << "): ";
        for (auto p : candidates_pos) std::cout << p << " ";
        std::cout << "\n";
        for (auto p : candidates_pos) {
            if (p + q.size() <= ref.size()) {
                std::cout << "  Candidate seq @ " << p << ": "
                          << ref.substr(p, q.size()) << "\n";
            }
        }

        // === metrics ===
        auto truth_str = posVecToStrVec(truth_pos);
        auto cand_str  = posVecToStrVec(candidates_pos);

        metrics::report(truth_str, cand_str);
        totalTP += metrics::countTP(truth_str, cand_str);
        totalFP += metrics::countFP(truth_str, cand_str);
        totalFN += metrics::countFN(truth_str, cand_str);
    }

    // ----- overall -----
    std::cout << "\n===== Overall Metrics =====\n";
    std::cout << "Total TP: " << totalTP << "\n";
    std::cout << "Total FP: " << totalFP << "\n";
    std::cout << "Total FN: " << totalFN << "\n";
    double recall = (totalTP + totalFN > 0) ? double(totalTP) / double(totalTP + totalFN) : 0.0;
    double precision = (totalTP + totalFP > 0) ? double(totalTP) / double(totalTP + totalFP) : 0.0;
    std::cout << "Overall Recall: " << recall << "\n";
    std::cout << "Overall Precision: " << precision << "\n";

    return 0;
}
