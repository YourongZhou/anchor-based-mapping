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
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/align.h>

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
std::pair<vector<int>, size_t> retrieveCandidates_mtree(
    MTree &mtree,
    const std::string &query,
    int maxDist)
{
    std::vector<int> results;

    Substring query_anchor{query, -1}; // pos = -1 表示查询
    
    // 1. 调用函数，'result_struct' 是包含 .matches 和 .nodeAccesses 的结构体
    auto result_struct = mtree.get_nearest_by_range(query_anchor, maxDist);

    // 2. 遍历结构体内部的 .matches 向量
    //    (您也可以使用更简洁的 C++11 范围for循环)
    for (const auto &res : result_struct.matches) {
        // 'res' 现在直接是 result_item (包含 .data 和 .distance)
        results.push_back(res.data.pos);
    }

    // 3. (可选) 获取节点访问次数
    size_t accesses = result_struct.nodeAccesses;
    std::cout << "M-Tree search node accesses: " << accesses << std::endl;

    return {results, accesses};
}


std::vector<int> retrieveCandidates_sae(
    seqan::Index<seqan::Dna5String, seqan::FMIndex<>> &fm_index, // 显式类型
    const seqan::Dna5String &ref_seq,                           // 显式类型
    const std::string &query,
    int maxDist,
    int seed_len) 
{
    // 显式使用 seqan::Score 和 seqan::Simple
    seqan::Score<int, seqan::Simple> scoring(1, -1, -1); // match=+1, mismatch=-1, gap=-1
    int xDropThreshold = maxDist * 4; // 经验值
    // ======================

    std::set<int> unique_positions;
    std::vector<int> results;

    std::string qseq_str = query;
    std::transform(qseq_str.begin(), qseq_str.end(), qseq_str.begin(), ::toupper);
    
    // 显式使用 seqan::Dna5String 和 seqan::assign
    seqan::Dna5String qseq;
    seqan::assign(qseq, qseq_str);

    // 显式类型定义
    typedef seqan::Seed<seqan::Simple> TSeed;
    
    // 显式使用 seqan::Finder
    seqan::Finder<seqan::Index<seqan::Dna5String, seqan::FMIndex<>>> finder(fm_index);

    // 遍历所有可能的种子
    for (size_t i = 0; i + seed_len <= qseq_str.size(); i += 3) {
        std::string seed_str = qseq_str.substr(i, seed_len);
        seqan::Dna5String seed;
        seqan::assign(seed, seed_str);
        
        seqan::clear(finder);
        
        // 显式使用 seqan::find, seqan::position
        while (seqan::find(finder, seed)) {
            size_t seed_ini_pos_on_ref = seqan::position(finder);
            
            // 显式使用 TSeed 构造函数
            TSeed s(seed_ini_pos_on_ref, i, 
                    seed_ini_pos_on_ref + seed_len - 1, i + seed_len - 1);
            
            // 显式使用 seqan::extendSeed, seqan::EXTEND_BOTH, seqan::GappedXDrop
            seqan::extendSeed(s, ref_seq, qseq, seqan::EXTEND_BOTH, scoring, xDropThreshold, seqan::GappedXDrop());

            // 显式使用 seqan::beginPositionH
            unsigned sb = seqan::beginPositionH(s);
            
            // 存储起始位置
            unique_positions.insert((int)sb); 
        }
    }

    results.reserve(unique_positions.size());
    results.assign(unique_positions.begin(), unique_positions.end());
    
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
    string fastq_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fastq/ERR15404863_1.fastq";
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
    int min_node_capacity = 50;              // min node capacity
    int max_node_capacity = -1;             // max node capacity (-1 表示默认/自动)
    double leaf_radius_threshold = 0.0;     // 叶子节点半径允许的阈值
    int compactness_min_capacity = 25;      // 触发紧凑性分裂的最小节点容量
    double compactness_radius_factor = 1.1; // 触发紧凑性分裂的半径膨胀比例
    std::string split_strategy_name = "mtree"; // 默认使用 TwoWaySplitStrategy
    bool auto_gen = false;
    int seed_len = 10;
    string method = "mtree"; // 方法，mtree, seed-and-extend, 还是 anchor

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
        } else if (arg == "--min_node_capacity" && i + 1 < argc) {
            min_node_capacity = stoi(argv[++i]);
        } else if (arg == "--max_node_capacity" && i + 1 < argc) {
            max_node_capacity = stoi(argv[++i]);
        } else if (arg == "--leaf_radius_threshold" && i + 1 < argc) {
            leaf_radius_threshold = stod(argv[++i]);
        } else if (arg == "--compactness_min_capacity" && i + 1 < argc) {
            compactness_min_capacity = stoi(argv[++i]);
        } else if (arg == "--compactness_radius_factor" && i + 1 < argc) {
            compactness_radius_factor = stod(argv[++i]);
        } else if (arg == "--split_strategy" && i + 1 < argc) {
            split_strategy_name = argv[++i];
            if (split_strategy_name != "mtree" && split_strategy_name != "fmtree") {
                cerr << "Invalid value for --split_strategy (expect mtree or fmtree)\n";
                return 1;
            }
        } else if (arg == "--auto_gen" && i + 1 < argc) {
            string val = argv[++i];
            if (val == "true" || val == "1") {
                auto_gen = true;
            } else if (val == "false" || val == "0") {
                auto_gen = false;
            } else {
                cerr << "Invalid value for --auto_gen (expect true/false or 1/0)\n";
                return 1;
            }
        } else if (arg == "--seed_len" && i + 1 < argc) {
            seed_len = stod(argv[++i]);
        } else if (arg == "--method" && i + 1 < argc) {
            method = argv[++i];
            if (method != "mtree" && method != "sae" && method != "anchor") {
                cerr << "Invalid value for --split_strategy (expect mtree or fmtree)\n";
                return 1;
            }
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
    cout << "--- M-Tree Config ---\n";
    cout << "min_node_capacity: " << min_node_capacity << "\n";
    cout << "max_node_capacity: " << max_node_capacity << "\n";
    cout << "leaf_radius_threshold: " << leaf_radius_threshold << "\n";
    cout << "compactness_min_capacity: " << compactness_min_capacity << "\n";
    cout << "compactness_radius_factor: " << compactness_radius_factor << "\n";
    cout << "auto generate queries: " << auto_gen << "\n";
    rng.seed(seed);
    time_t time_start, time_read, time_tree, time_overlap, time_truth, time_query;
    time_start = time(NULL);
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
    time_read = time(NULL);

    // ----- 生成 m-tree -----
    MTree mtree(
        min_node_capacity,          // min node capacity
        max_node_capacity,          // max node capacity
        leaf_radius_threshold,      // 叶子节点半径允许的阈值
        compactness_min_capacity,   // 触发紧凑性分裂的最小节点容量
        compactness_radius_factor,  // 触发紧凑性分裂的半径膨胀比例
        Distance(),                 // 距离函数
        SplitStrategyType()         // split function
    );

    if (method == "mtree"){
        build_mtree_from_ref(ref, anchor_len, mtree);
        time_tree = time(NULL);
        // 输出 mtree overlap
        // mtree.print_overlap_info();

        // 检查树状态
        cout << "Validating tree structure:" << endl;
        // mtree._check();
        cout << "Validation successful." << endl;
        time_overlap = time(NULL);   // 输出mtree各层次半径
        // print_mtree_radius_distribution(mtree);
    }

    // ----- 生成 anchors -----
    auto anchors = generate_anchors_from_seq(ref, "ref1", anchor_len, num_anchors);
    cout << "Generated " << anchors.size() << " anchors\n";
    // 用 set 检查去重后数量
    unordered_set<string> uniq;
    for (auto &a : anchors) {
        uniq.insert(a.seq);   // 注意这里用 seq，而不是 id
    }

    cout << "Unique anchors: " << uniq.size() << endl;

    // 如果有重复，提示
    if (uniq.size() < anchors.size()) {
        cout << "Warning: Found " 
                  << (anchors.size() - uniq.size()) 
                  << " duplicate anchors!" << endl;
    }

    unordered_map<string, vector<pair<int,int>>> anchor_index;
    if (method == "anchor"){
        // ----- 构建 anchor index -----
        cout << "Building anchor index:\n";
        auto anchor_index = build_anchor_index(anchors, ref, anchor_len);
    }

    // ----- 生成 queries -----
    vector<string> queries;
    if (auto_gen){
        queries = simulate_queries(ref, num_queries, query_len);
        cout << "Simulated " << num_queries << " queries\n";
    } else {
        try {
                vector<FastqRecord> fastq_records = read_fastq_seqan(fastq_path);
                
                int count = 0;
                for (const auto& record : fastq_records) {
                    if (count++ >= num_queries){
                        break;
                    }
                    queries.push_back(record.seq.substr(0, query_len));
                }
                cout << "Read " << queries.size() << " queries from FASTQ file: " << fastq_path << "\n";
            } catch (const runtime_error& e) {
                cerr << "Error reading FASTQ: " << e.what() << "\n";
            }
    }

    // // 只保留周围某个距离内有 anchor 的 query
    // queries = filter_queries_by_anchor_index(queries, anchor_index, 3);
    // num_queries = queries.size();

    // ----- ground truth -----
    cout << "Generating ground truths:\n";
    cout << "Number of queries: " << num_queries << endl;
    vector<vector<int>> truth_positions;
    truth_positions.reserve(num_queries);
    for (const auto &q : queries) {
        truth_positions.push_back(find_all_occurrences_approx(q, ref, maxDist));
    }
    time_truth = time(NULL);
    
    seqan::Index<seqan::Dna5String, seqan::FMIndex<>> fm_index;
    seqan::Dna5String ref_seq;
    if (method == "sae"){
        // ----- FMindex -----
        seqan::assign(ref_seq, ref);
        // 使用 Dna5String 作为文本类型，FMIndex<> 作为索引类型
        seqan::Index<seqan::Dna5String, seqan::FMIndex<>> fm_index(ref_seq);
        try {
            seqan::indexRequire(fm_index, seqan::FibreSA());
            std::cout << "[INFO] FM-index constructed successfully and Suffix Array loaded.\n";
        } catch (const seqan::Exception &e) {
            std::cerr << "[ERROR] Failed to construct or load Suffix Array for FM-index: " << e.what() << "\n";
            return 1; // 或者采取其他错误处理措施
        }
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

        
        // === metrics ===
        size_t count = 0;  // 保存 dist_list.size()
        const auto &truth_str = posVecToStrVec(truth_positions[i]);
        cout << "truth positions: " << endl;
        for (size_t j = 0; j < truth_positions[i].size(); ++j) {
            std::cout << j << ": " << truth_positions[i][j] << std::endl;
        }

        vector<string> cand_str;
        if (method == "mtree"){
            auto [pos, access] = retrieveCandidates_mtree(mtree, queries[i], maxDist);
            cand_str = posVecToStrVec(pos);
        } else if (method == "sae"){
            cand_str = posVecToStrVec(retrieveCandidates_sae(fm_index, ref_seq, queries[i], maxDist, seed_len));
        } else{
            cand_str = posVecToStrVec(retrieveCandidates_anchor(queries[i], anchor_index, anchors[0].seq.size(), maxDist, use_all, start_idx, end_idx, last_random, anchor_radius = anchor_radius, &count));
        };

        dist_counts[i] = count;
        cout << "candidate positions: " << endl;
        for (size_t j = 0; j < cand_str.size(); ++j) {
            std::cout << j << ": " << cand_str[j] << std::endl;
        }

        int tp = metrics::countTP(truth_str, cand_str);
        int fp = metrics::countFP(truth_str, cand_str);
        int fn = metrics::countFN(truth_str, cand_str);
        auto [avg_dist, max_dist] = metrics::evaluateDistances(queries[i], cand_str, ref);

        // double recall = (tp + fn > 0) ? double(tp) / double(tp + fn) : 0.0;
        // double precision = (tp + fp > 0) ? double(tp) / double(tp + fp) : 0.0;
        // double fp_over_tp = (tp > 0) ? double(fp) / double(tp) : 0.0;

        sumTP += tp;
        sumFP += fp;
        sumFN += fn;
        // sumRecall += recall;
        // sumPrecision += precision;
        // sum_fp_over_tp += fp_over_tp;
        sum_avg_dist += avg_dist;
        sum_max_dist += max_dist;
        // metrics::report(truth_str, cand_str);
    }
    time_query = time(NULL);
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

    // 总指标
    double recall = (sumTP + sumFN > 0) ? double(sumTP) / double(sumTP + sumFN) : 0.0;
    double precision = (sumTP + sumFP > 0) ? double(sumTP) / double(sumTP + sumFP) : 0.0;
    double fp_over_tp = (sumTP > 0) ? double(sumFP) / double(sumTP) : 0.0;
    double avg_avg_dist = sum_avg_dist / num_queries;
    double avg_max_dist = sum_max_dist / num_queries;

    cout << "\n===== Average per-query Metrics =====\n";
    cout << "Average TP: " << (double)sumTP / num_queries << "\n";
    cout << "Average FP: " << (double)sumFP / num_queries << "\n";
    cout << "Average FN: " << (double)sumFN / num_queries << "\n";
    cout << "Average Recall: " << recall << "\n";
    cout << "Average Precision: " << precision << "\n";
    cout << "Average FP/TP: " << fp_over_tp << "\n";
    cout << "Average average distance: " << avg_avg_dist << "\n";
    cout << "Average maximum distance: " << avg_max_dist << "\n";

    cout << "===== Average time lapse =====";
    printf("\nRead FASTQA file time:%ld", (time_read-time_start));
    printf("\nBuild Tree time:%ld", (time_tree-time_read));
    printf("\nCalculate overlap time:%ld", (time_overlap - time_tree));
    printf("\nGet ground truth time:%ld", (time_truth - time_overlap));
    printf("\nQuery time:%ld", (time_query - time_truth));

    return 0;
}
