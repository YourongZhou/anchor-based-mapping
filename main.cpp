// main.cpp

#include "fasta_utils_seqan.hpp"
#include "fastq_utils_seqan.hpp"
#include "levenshtein.hpp"
#include "metrics.h"

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

// 兼容：你之前有一个 retrieveCandidates(anchorDist, maxDist)；
// 这里我们实现一个新的签名：给定 query, anchors, ref -> 返回 candidate positions（int vector）
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

// 小工具：把 int positions 转成 string id （用于复用现有 string-based metrics）
inline std::vector<std::string> posVecToStrVec(const std::vector<int> &pos) {
    std::vector<std::string> out;
    out.reserve(pos.size());
    for (int p : pos) out.push_back(std::to_string(p));
    return out;
}

// -----------------------------
// main demo
// -----------------------------
int main() {
    // ----- 参数（你可以按需修改） -----
    const std::string fasta_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fasta/ecoli.fa";
    const size_t truncate_ref_len = 10000;   // 只取前 10000 个碱基
    const size_t anchor_len = 20;            // anchor 长度为 20
    const int num_anchors = 50000;             // anchor 个数
    const int num_queries = 100;             // 模拟 queries 数量
    const size_t query_len = 20;             // query 长度（同 anchor_len）
    const int maxDist = 2;                   // anchor 与 query 的最大距离阈值

    // ----- 读取 FASTA 并截断到前 10000bp -----
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

    // ----- 生成 anchors（从单条 reference 里随机采样） -----
    auto anchors = generate_anchors_from_seq(ref, "ref1", anchor_len, num_anchors);
    std::cout << "Generated " << anchors.size() << " anchors\n";

    // ----- 生成模拟 queries（从 ref 随机抽取） -----
    auto queries = simulate_queries(ref, num_queries, query_len);
    std::cout << "Simulated " << queries.size() << " queries\n";

    // ----- 计算 ground truth（每个 query 在 ref 上的所有位置） -----
    std::vector<std::vector<int>> truth_positions;
    truth_positions.reserve(queries.size());
    for (const auto &q : queries) {
        auto pos = find_all_occurrences(q, ref);
        // pos 可能为空（如果 query 没在 ref 出现，这里理论上不会发生因为 queries 源自 ref）
        truth_positions.push_back(std::move(pos));
    }

    // ----- 对每个 query 用 anchor-based 方法检索 candidate positions，并评估 -----
        // ----- 对每个 query 用 anchor-based 方法检索 candidate positions，并评估 -----
    int totalTP = 0, totalFP = 0, totalFN = 0;
    for (size_t i = 0; i < queries.size(); ++i) {
        const auto &q = queries[i];
        auto truth_pos = truth_positions[i];                 // vector<int>
        auto candidates_pos = retrieveCandidates(q, anchors, ref, maxDist); // vector<int>

        // 为了复用现有 metrics（接受 vector<string>），把位置转换为 string
        auto truth_str = posVecToStrVec(truth_pos);
        auto cand_str  = posVecToStrVec(candidates_pos);

        // ====== 打印详细信息 ======
        std::cout << "\n=== Query [" << i << "] ===\n";
        std::cout << "Query seq: " << q << "\n";

        std::cout << "Truth positions (" << truth_pos.size() << "): ";
        for (auto p : truth_pos) std::cout << p << " ";
        std::cout << "\n";
        // 打印 truth 的序列
        for (auto p : truth_pos) {
            if (p + q.size() <= ref.size()) {
                std::cout << "  Truth seq @ " << p << ": "
                          << ref.substr(p, q.size()) << "\n";
            }
        }

        std::cout << "Candidate positions (" << candidates_pos.size() << "): ";
        for (auto p : candidates_pos) std::cout << p << " ";
        std::cout << "\n";
        // 打印 candidate 的序列
        for (auto p : candidates_pos) {
            if (p + q.size() <= ref.size()) {
                std::cout << "  Candidate seq @ " << p << ": "
                          << ref.substr(p, q.size()) << "\n";
            }
        }

        // ====== metrics ======
        metrics::report(truth_str, cand_str);
        totalTP += metrics::countTP(truth_str, cand_str);
        totalFP += metrics::countFP(truth_str, cand_str);
        totalFN += metrics::countFN(truth_str, cand_str);
    }


    // ----- overall metrics -----
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
