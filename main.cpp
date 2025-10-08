// main.cpp

#include "rng.h"
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
#include <cstdlib>
#include <nlohmann/json.hpp>


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

    // static std::mt19937 rng(std::random_device{}());
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
    int max_dist_query, // 扩展窗口大小（d +/- max_dist_query）
    bool use_all = true, // 是否选择全部 anchor
    int start_idx = 0,  // 起始比例
    int end_idx   = 20,  // 结束比例
    bool last_random = false,  // 最后一个是否随机选择（当按顺序选择的时候） 
    int anchor_radius = 3
) {
    int qlen = (int)query.size();
    if (qlen < k) return {};

    // ---------- Step 1: 计算 query 与所有 anchor 的距离 ----------
    std::vector<std::pair<std::string,int>> dist_list;
    dist_list.reserve(anchor_index.size());

    for (const auto &kv : anchor_index) {
        const std::string &anchor_seq = kv.first;
        if ((int)anchor_seq.size() != k) continue;
        int d = min_edit_distance_window(anchor_seq, query);

        // ✅ 只保留距离在 anchor_radius 以内的 anchor
        if (d <= anchor_radius)
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

    std::vector<std::unordered_set<int>> candidate_sets;
    candidate_sets.reserve(64);

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
    std::cout << "\nAverage anchor - query distance: " << anchor_query_dist / selected.size() << " " << selected.size();

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

int main(int argc, char* argv[]) {
    // ----- 参数 -----
    std::string fasta_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fasta/ecoli.fa";
    size_t truncate_ref_len = 10000;    // 只取前 10000 个碱基
    size_t anchor_len = 20;             // anchor 长度
    int num_anchors = 100;              // anchor 个数
    int num_queries = 100;              // query 数量
    size_t query_len = anchor_len;      // query 长度
    int maxDist = 3;                    // 容差距离
    bool use_all = true; // 是否选择全部 anchor
    double start_idx = 0.0;  // 起始比例
    double end_idx   = 0.2;  // 结束比例
    bool last_random = false; // 最后一个随机选
    unsigned int seed = 42; // 默认 seed
    int anchor_radius = 3; // 选择的 anchor 到 query 的距离

    // ----- 解析命令行参数 -----
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--truncate_ref_len" && i + 1 < argc) {
            truncate_ref_len = std::stoul(argv[++i]);
        } else if (arg == "--anchor_len" && i + 1 < argc) {
            anchor_len = std::stoul(argv[++i]);
        } else if (arg == "--num_anchors" && i + 1 < argc) {
            num_anchors = std::stoi(argv[++i]);
        } else if (arg == "--num_queries" && i + 1 < argc) {
            num_queries = std::stoi(argv[++i]);
        } else if (arg == "--maxDist" && i + 1 < argc) {
            maxDist = std::stoi(argv[++i]);
        } else if (arg == "--use_all" && i + 1 < argc) {
            std::string val = argv[++i];
            if (val == "true" || val == "1") {
                use_all = true;
            } else if (val == "false" || val == "0") {
                use_all = false;
            } else {
                std::cerr << "Invalid value for --use_all (expect true/false or 1/0)\n";
                return 1;
            }
        } else if (arg == "--start_idx" && i + 1 < argc) {
            start_idx = std::stod(argv[++i]);
        } else if (arg == "--end_idx" && i + 1 < argc) {
            end_idx = std::stod(argv[++i]);
        } else if (arg == "--last_random" && i + 1 < argc) {
            std::string val = argv[++i];
            if (val == "true" || val == "1") {
                last_random = true;
            } else if (val == "false" || val == "0") {
                last_random = false;
            } else {
                std::cerr << "Invalid value for --last_random (expect true/false or 1/0)\n";
                return 1;
            }
        } else if (arg == "--seed" && i + 1 < argc) {
            seed = std::stoul(argv[++i]);
        } else if (arg == "--anchor_radius" && i + 1 < argc) {
            end_idx = std::stod(argv[++i]);
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            return 1;
        }
    }

    // ----- 打印参数确认 -----
    std::cout << "fasta_path: " << fasta_path << "\n";
    std::cout << "truncate_ref_len: " << truncate_ref_len << "\n";
    std::cout << "anchor_len: " << anchor_len << "\n";
    std::cout << "num_anchors: " << num_anchors << "\n";
    std::cout << "num_queries: " << num_queries << "\n";
    std::cout << "maxDist: " << maxDist << "\n";
    std::cout << "use_all: " << (use_all ? "true" : "false") << "\n";
    std::cout << "start_idx: " << start_idx << "\n";
    std::cout << "end_idx: " << end_idx << "\n";
    std::cout << "last_random: " << last_random << "\n";
    std::cout << "anchor radius: " << anchor_radius << "\n"; 
    std::cout << "seed: " << seed << "\n"; 
    rng.seed(seed);

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
    // 用 set 检查去重后数量
    std::unordered_set<std::string> uniq;
    for (auto &a : anchors) {
        uniq.insert(a.seq);   // 注意这里用 seq，而不是 id
    }

    std::cout << "Unique anchors: " << uniq.size() << std::endl;

    // 如果有重复，提示
    if (uniq.size() < anchors.size()) {
        std::cout << "Warning: Found " 
                  << (anchors.size() - uniq.size()) 
                  << " duplicate anchors!" << std::endl;
    }

    // ----- 生成 queries -----
    auto queries = simulate_queries(ref, num_queries, query_len);
    std::cout << "Simulated " << queries.size() << " queries\n";

    // ----- ground truth -----
    std::cout << "Generating ground truths:\n";
    std::vector<std::vector<int>> truth_positions;
    truth_positions.reserve(queries.size());
    for (const auto &q : queries) {
        truth_positions.push_back(find_all_occurrences_approx(q, ref, maxDist));
    }

    // ----- 构建 anchor index -----
    std::cout << "Building anchor index:\n";
    auto anchor_index = build_anchor_index(anchors, ref, anchor_len);
    using json = nlohmann::json;

    json j;
    for (auto &[seq, matches] : anchor_index) {
        for (auto &[pos, dist] : matches) {
            j["anchors"].push_back({{"seq", seq}, {"pos", pos}, {"dist", dist}});
        }
    }
    j["queries"] = queries;

    std::ofstream out("anchor_index.json");
    out << j.dump(2);


    // ----- 查询 + 评估 -----
    // 用于累加每条 query 的 recall / precision
    int sumTP = 0;
    int sumFP = 0;
    int sumFN = 0;
    double sumRecall = 0.0;
    double sumPrecision = 0.0;
    double sum_fp_over_tp = 0.0;
    double sum_avg_dist = 0.0;
    double sum_max_dist = 0.0;
    for (size_t i = 0; i < queries.size(); ++i) {
        const auto &q = queries[i];
        const auto &truth_pos = truth_positions[i];

        // === 用 anchor-based 检索 ===
        // auto candidates_pos = retrieveCandidates_anchor(q, anchor_index, anchors[0].seq.size(), maxDist, use_all, start_idx, end_idx);

        // // === 打印 ===
        // std::cout << "\n=== Query [" << i << "] ===\n";
        // std::cout << "Query seq: " << q << "\n";

        // std::cout << "Truth positions (" << truth_pos.size() << "): ";
        // for (auto p : truth_pos) std::cout << p << " ";
        // std::cout << "\n";
        // for (auto p : truth_pos) {
        //     if (p + q.size() <= ref.size()) {
        //         std::cout << "  Truth seq @ " << p << ": "
        //                   << ref.substr(p, q.size()) << "\n";
        //     }
        // }

        // std::cout << "Candidate positions (" << candidates_pos.size() << "): ";
        // for (auto p : candidates_pos) std::cout << p << " ";
        // std::cout << "\n";
        // for (auto p : candidates_pos) {
        //     if (p + q.size() <= ref.size()) {
        //         std::cout << "  Candidate seq @ " << p << ": "
        //                   << ref.substr(p, q.size()) << "\n";
        //     }
        // }

        // === metrics ===
        const auto &truth_str = posVecToStrVec(truth_positions[i]);
        const auto &cand_str  = posVecToStrVec(retrieveCandidates_anchor(queries[i], anchor_index, anchors[0].seq.size(), maxDist, use_all, start_idx, end_idx, last_random));

        int tp = metrics::countTP(truth_str, cand_str);
        int fp = metrics::countFP(truth_str, cand_str);
        int fn = metrics::countFN(truth_str, cand_str);
        auto [avg_dist, max_dist] = metrics::evaluateDistances(queries[i], cand_str);


        double recall = (tp + fn > 0) ? double(tp) / double(tp + fn) : 0.0;
        double precision = (tp + fp > 0) ? double(tp) / double(tp + fp) : 0.0;
        double fp_over_tp = double(fp) / double(tp);


        sumTP += tp;
        sumFP += fp;
        sumFN += fn;
        sumRecall += recall;
        sumPrecision += precision;
        sum_fp_over_tp += fp_over_tp;
        sum_avg_dist += avg_dist;
        sum_max_dist += max_dist;
        // metrics::report(truth_str, cand_str);
    }

    // ----- overall -----
    std::cout << "\n===== Experiment Parameters =====\n";
    std::cout << "Reference length truncated to: " << truncate_ref_len << " bp\n";
    std::cout << "Anchor/Query length: " << anchor_len << "\n";
    std::cout << "Number of anchors: " << num_anchors << "\n";
    std::cout << "Number of queries: " << num_queries << "\n";
    std::cout << "Maximum allowed distance: " << maxDist << "\n";

    // 平均每条 query 的指标
    double avgRecall = sumRecall / queries.size();
    double avgPrecision = sumPrecision / queries.size();

    std::cout << "\n===== Average per-query Metrics =====\n";
    std::cout << "Average TP: " << (double)sumTP / queries.size() << "\n";
    std::cout << "Average FP: " << (double)sumFP / queries.size() << "\n";
    std::cout << "Average FN: " << (double)sumFN / queries.size() << "\n";
    std::cout << "Average Recall: " << avgRecall << "\n";
    std::cout << "Average Precision: " << avgPrecision << "\n";
    std::cout << "Average FP/TP: " << sum_fp_over_tp << "\n";
    std::cout << "Average average distance: " << sum_avg_dist << "\n";
    std::cout << "Average maximum distance: " << sum_max_dist << "\n";

    return 0;
}
