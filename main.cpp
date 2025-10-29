// main.cpp

#include "rng.h"
#include "fasta_utils_seqan.hpp"
#include "fastq_utils_seqan.hpp"
#include "levenshtein.hpp"
#include "metrics.h"
#include "anchor_gen.h"
#include "tools.h"
#include "candidate_retriever.h"
#include "index_query.h"
#include "mtree.h"
#include "functions.h"

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

using namespace std;
using namespace mt;


// Levenshtein 距离封装，只比较 seq
struct SubstringLevDist {
    double operator()(const Substring &a, const Substring &b) const {
        return static_cast<double>(levenshtein(a.seq, b.seq));
    }
};

using Data = Substring;
using Distance = SubstringLevDist;
using CachedDistance = functions::cached_distance_function<Data, Distance>;

// split function 类型定义
// using SplitStrategyType = functions::TwoWaySplitStrategy<
//     functions::random_promotion,
//     functions::balanced_partition
// >;
using SplitStrategyType = functions::OptimizedKSplitStrategy;   
// using SplitFunc = functions::split_function<
//     functions::random_promotion,
//     functions::balanced_partition>;

// 定义最终 MTree 类型
using MTree = mtree<Data, Distance, SplitStrategyType>;

// ======================== 构建函数 ========================
void build_mtree_from_ref(const std::string &ref, int k, MTree &tree) {
    std::unordered_set<std::string> inserted;  // 记录已插入的序列
    size_t count = 0;

    for (size_t i = 0; i + k <= ref.size(); ++i) {
        std::string subseq = ref.substr(i, k);
        if (inserted.find(subseq) != inserted.end())
            continue; // 已插入则跳过

        Substring a{subseq, static_cast<int>(i)};
        tree.add(a);
        inserted.insert(subseq);

        count++;
        if (count % 1000 == 0)
            std::cout << "Inserted " << count << " unique anchors into MTree." << std::endl;
    }

    std::cout << "MTree construction completed. Total unique anchors: " << count << std::endl;
}

// ======================== 查询函数 ========================
vector<int> retrieveCandidates_mtree(
    MTree &mtree,
    const std::string &query,
    int maxDist)
{
    std::vector<int> results;

    Substring query_anchor{query, -1}; // pos = -1 表示查询
    auto matches = mtree.get_nearest_by_range(query_anchor, maxDist);

    for (auto it = matches.begin(); it != matches.end(); ++it) {
        const auto &res = *it;  // res 包含 (obj, distance)
        results.push_back(res.data.pos);
    }
    return results;
}

// // 递归遍历所有节点
// template <typename Data, typename Distance, typename SplitFunc>
// void traverse_mtree_levels(
//     typename mtree<Data, Distance, SplitFunc>::Node* node,
//     int level,
//     std::map<int, std::vector<float>>& radii_by_level
// ) {
//     if (!node) return;

//     // 遍历该节点的每个 entry
//     for (const auto& entry : node->entries) {
//         // 记录该 entry 的覆盖半径
//         radii_by_level[level].push_back(entry.radius);

//         // 若该 entry 有子节点，则递归下去
//         if (entry.child != nullptr) {
//             traverse_mtree_levels<Data, Distance, SplitFunc>(entry.child, level + 1, radii_by_level);
//         }
//     }
// }

// // 打印统计信息
// template <typename Data, typename Distance, typename SplitFunc>
// void print_mtree_radius_distribution(const mtree<Data, Distance, SplitFunc>& tree) {
//     std::map<int, std::vector<float>> radii_by_level;
//     traverse_mtree_levels<Data, Distance, SplitFunc>(tree.root, 0, radii_by_level);

//     std::cout << "\n===== M-tree Radius Distribution =====" << std::endl;
//     for (const auto& [level, radii] : radii_by_level) {
//         double mean = 0, maxr = 0, minr = 1e9;
//         for (auto r : radii) {
//             mean += r;
//             maxr = std::max(maxr, (double)r);
//             minr = std::min(minr, (double)r);
//         }
//         mean /= radii.size();
//         std::cout << "Level " << level
//                   << " | Nodes: " << radii.size()
//                   << " | mean radius = " << mean
//                   << " | min = " << minr
//                   << " | max = " << maxr
//                   << std::endl;
//     }
// }


int main(int argc, char* argv[]) {
    // ----- 参数 -----
    string fasta_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fasta/ecoli.fa";
    size_t truncate_ref_len = 10000;    // 只取前 10000 个碱基
    size_t anchor_len = 20;             // anchor 长度
    int num_anchors = 100;              // anchor 个数
    int num_queries = 100;              // query 数量
    int maxDist = 3;                    // 容差距离
    bool use_all = true; // 是否选择全部 anchor
    double start_idx = 0.0;  // 起始比例
    double end_idx   = 0.2;  // 结束比例
    bool last_random = false; // 最后一个随机选
    unsigned int seed = 42; //  默认 seed
    bool use_anchor_radius = false; // 默认不用某个距离内的 anchor
    int anchor_radius = 3; // 选择的 anchor 到 query 的距离

    // ----- 解析命令行参数 -----
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--truncate_ref_len" && i + 1 < argc) {
            truncate_ref_len = stoul(argv[++i]);
        } else if (arg == "--anchor_len" && i + 1 < argc) {
            anchor_len = stoul(argv[++i]);
        } else if (arg == "--num_anchors" && i + 1 < argc) {
            num_anchors = stoi(argv[++i]);
        } else if (arg == "--num_queries" && i + 1 < argc) {
            num_queries = stoi(argv[++i]);
        } else if (arg == "--maxDist" && i + 1 < argc) {
            maxDist = stoi(argv[++i]);
        } else if (arg == "--use_all" && i + 1 < argc) {
            string val = argv[++i];
            if (val == "true" || val == "1") {
                use_all = true;
            } else if (val == "false" || val == "0") {
                use_all = false;
            } else {
                cerr << "Invalid value for --use_all (expect true/false or 1/0)\n";
                return 1;
            }
        } else if (arg == "--start_idx" && i + 1 < argc) {
            start_idx = stod(argv[++i]);
        } else if (arg == "--end_idx" && i + 1 < argc) {
            end_idx = stod(argv[++i]);
        } else if (arg == "--last_random" && i + 1 < argc) {
            string val = argv[++i];
            if (val == "true" || val == "1") {
                last_random = true;
            } else if (val == "false" || val == "0") {
                last_random = false;
            } else {
                cerr << "Invalid value for --last_random (expect true/false or 1/0)\n";
                return 1;
            }
        } else if (arg == "--seed" && i + 1 < argc) {
            seed = stoul(argv[++i]);
        } else if (arg == "--anchor_radius" && i + 1 < argc) {
            anchor_radius = stod(argv[++i]);
        } else {
            cerr << "Unknown or incomplete argument: " << arg << "\n";
            return 1;
        }
    }
    size_t query_len = anchor_len;      // query 长度

    // ----- 打印参数确认 -----
    cout << "fasta_path: " << fasta_path << "\n";
    cout << "truncate_ref_len: " << truncate_ref_len << "\n";
    cout << "anchor_len: " << anchor_len << "\n";
    cout << "num_anchors: " << num_anchors << "\n";
    cout << "num_queries: " << num_queries << "\n";
    cout << "maxDist: " << maxDist << "\n";
    cout << "use_all: " << (use_all ? "true" : "false") << "\n";
    cout << "start_idx: " << start_idx << "\n";
    cout << "end_idx: " << end_idx << "\n";
    cout << "last_random: " << last_random << "\n";
    cout << "anchor radius: " << anchor_radius << "\n"; 
    cout << "seed: " << seed << "\n"; 
    rng.seed(seed);

    // ----- 读取 FASTA 并截断 -----
    vector<FastaRecord> records;
    try {
        records = read_fasta_seqan(fasta_path);
    } catch (runtime_error &e) {
        cerr << "Error reading FASTA: " << e.what() << "\n";
        return 1;
    }
    if (records.empty()) {
        cerr << "No sequences in FASTA\n";
        return 1;
    }
    string ref = records[0].seq.substr(0, min(records[0].seq.size(), truncate_ref_len));
    cout << "Using reference (truncated) length = " << ref.size() << "\n";

    // ----- 生成 m-tree -----
    MTree mtree(
        2,              // min node capacity
        -1,             // max node capacity
        5,
        5,
        0.2,
        Distance(),     // 距离函数
        SplitStrategyType()  // split function
    );

    build_mtree_from_ref(ref, anchor_len, mtree);
    // 输出mtree各层次半径
    // print_mtree_radius_distribution(mtree);

    // // ----- 生成 anchors -----
    // auto anchors = generate_anchors_from_seq(ref, "ref1", anchor_len, num_anchors);
    // cout << "Generated " << anchors.size() << " anchors\n";
    // // 用 set 检查去重后数量
    // unordered_set<string> uniq;
    // for (auto &a : anchors) {
    //     uniq.insert(a.seq);   // 注意这里用 seq，而不是 id
    // }

    // cout << "Unique anchors: " << uniq.size() << endl;

    // // 如果有重复，提示
    // if (uniq.size() < anchors.size()) {
    //     cout << "Warning: Found " 
    //               << (anchors.size() - uniq.size()) 
    //               << " duplicate anchors!" << endl;
    // }

    // // ----- 构建 anchor index -----
    // cout << "Building anchor index:\n";
    // auto anchor_index = build_anchor_index(anchors, ref, anchor_len);
    // ----- 生成 queries -----
    auto queries = simulate_queries(ref, num_queries, query_len);
    cout << "Simulated " << num_queries << " queries\n";

    // // 只保留周围某个距离内有 anchor 的 query
    // queries = filter_queries_by_anchor_index(queries, anchor_index, 3);
    // num_queries = queries.size();

    // ----- ground truth -----
    cout << "Generating ground truths:\n";
    vector<vector<int>> truth_positions;
    truth_positions.reserve(num_queries);
    for (const auto &q : queries) {
        truth_positions.push_back(find_all_occurrences_approx(q, ref, maxDist));
    }

    // using json = nlohmann::json;

    // json j;
    // for (auto &[seq, matches] : anchor_index) {
    //     for (auto &[pos, dist] : matches) {
    //         j["anchors"].push_back({{"seq", seq}, {"pos", pos}, {"dist", dist}});
    //     }
    // }
    // j["queries"] = queries;

    // ofstream out("anchor_index.json");
    // out << j.dump(2);


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
    vector<size_t> dist_counts(num_queries);

    for (size_t i = 0; i < num_queries; ++i) {
        const auto &q = queries[i];
        const auto &truth_pos = truth_positions[i];

        // === 用 anchor-based 检索 ===
        // auto candidates_pos = retrieveCandidates_anchor(q, anchor_index, anchors[0].seq.size(), maxDist, use_all, start_idx, end_idx);

        // // === 打印 ===
        // cout << "\n=== Query [" << i << "] ===\n";
        // cout << "Query seq: " << q << "\n";

        // cout << "Truth positions (" << truth_pos.size() << "): ";
        // for (auto p : truth_pos) cout << p << " ";
        // cout << "\n";
        // for (auto p : truth_pos) {
        //     if (p + q.size() <= ref.size()) {
        //         cout << "  Truth seq @ " << p << ": "
        //                   << ref.substr(p, q.size()) << "\n";
        //     }
        // }

        // cout << "Candidate positions (" << candidates_pos.size() << "): ";
        // for (auto p : candidates_pos) cout << p << " ";
        // cout << "\n";
        // for (auto p : candidates_pos) {
        //     if (p + q.size() <= ref.size()) {
        //         cout << "  Candidate seq @ " << p << ": "
        //                   << ref.substr(p, q.size()) << "\n";
        //     }
        // }

        // === metrics ===
        size_t count = 0;  // 保存 dist_list.size()
        const auto &truth_str = posVecToStrVec(truth_positions[i]);
        const auto &cand_str  = posVecToStrVec(retrieveCandidates_mtree(mtree, queries[i], maxDist));
        // const auto &cand_str  = posVecToStrVec(retrieveCandidates_anchor(queries[i], anchor_index, anchors[0].seq.size(), maxDist, use_all, start_idx, end_idx, last_random, anchor_radius = anchor_radius, &count));
        dist_counts[i] = count;

        int tp = metrics::countTP(truth_str, cand_str);
        int fp = metrics::countFP(truth_str, cand_str);
        int fn = metrics::countFN(truth_str, cand_str);
        auto [avg_dist, max_dist] = metrics::evaluateDistances(queries[i], cand_str, ref);

        double recall = (tp + fn > 0) ? double(tp) / double(tp + fn) : 0.0;
        double precision = (tp + fp > 0) ? double(tp) / double(tp + fp) : 0.0;
        double fp_over_tp = (tp > 0) ? double(fp) / double(tp) : 0.0;

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

    // 保存到文件
    ofstream ofs("dist_counts_radius_" + to_string(anchor_radius) + ".txt");

    for (auto c : dist_counts) {
        ofs << c << "\n";
    }
    ofs.close();

    // ----- overall -----
    cout << "\n===== Experiment Parameters =====\n";
    cout << "Reference length truncated to: " << truncate_ref_len << " bp\n";
    cout << "Anchor/Query length: " << anchor_len << "\n";
    cout << "Number of anchors: " << num_anchors << "\n";
    cout << "Number of queries: " << num_queries << "\n";
    cout << "Maximum allowed distance: " << maxDist << "\n";

    // 平均每条 query 的指标
    double avgRecall = sumRecall / num_queries;
    double avgPrecision = sumPrecision / num_queries;
    double avg_fp_over_tp = sum_fp_over_tp / num_queries;
    double avg_avg_dist = sum_avg_dist / num_queries;
    double avg_max_dist = sum_max_dist / num_queries;

    cout << "\n===== Average per-query Metrics =====\n";
    cout << "Average TP: " << (double)sumTP / num_queries << "\n";
    cout << "Average FP: " << (double)sumFP / num_queries << "\n";
    cout << "Average FN: " << (double)sumFN / num_queries << "\n";
    cout << "Average Recall: " << avgRecall << "\n";
    cout << "Average Precision: " << avgPrecision << "\n";
    cout << "Average FP/TP: " << avg_fp_over_tp << "\n";
    cout << "Average average distance: " << avg_avg_dist << "\n";
    cout << "Average maximum distance: " << avg_max_dist << "\n";

    return 0;
}
