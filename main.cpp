// main.cpp

#include "tools.h"
#include <omp.h>

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
    std::string index_dir = "./index_cache"; // 默认缓存目录
    bool force_rebuild = false;

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
        } else if (arg == "--index_dir" && i + 1 < argc) {
        index_dir = argv[++i];
        } else if (arg == "--force_rebuild") {
            force_rebuild = true;
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

    unordered_map<string, vector<pair<int,int>>> anchor_index;
    vector<FastaRecord> anchors;
    if (method == "anchor"){
        // ----- 生成 anchors -----
        anchors = generate_anchors_from_seq(ref, "ref1", anchor_len, num_anchors);
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
        // ----- 构建 anchor index -----
        cout << "Building anchor index:\n";
        anchor_index = build_anchor_index(anchors, ref, anchor_len);
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
    vector<vector<int>> truth_positions(num_queries); // 预先分配大小

    // 2. 并行化 Ground Truth 循环
    #pragma omp parallel for default(none) shared(queries, ref, maxDist, truth_positions, num_queries)
    for (size_t i = 0; i < num_queries; ++i) {
        truth_positions[i] = find_all_occurrences_approx(queries[i], ref, maxDist);
    }

    vector<string> filtered_queries;
    vector<vector<int>> filtered_truth_positions;
    size_t new_num_queries = 0;

    for (size_t i = 0; i < num_queries; ++i) {
        if (!truth_positions[i].empty()) {
            filtered_queries.push_back(queries[i]);
            filtered_truth_positions.push_back(truth_positions[i]);
            new_num_queries++;
        }
    }

    // 替换原始数据
    if (new_num_queries < num_queries) {
        std::cout << "\n[FILTER] Removed " << (num_queries - new_num_queries) 
                  << " queries with no ground truth (maxDist=" << maxDist << ").\n";
        queries = std::move(filtered_queries);
        truth_positions = std::move(filtered_truth_positions);
        num_queries = new_num_queries;
    } else {
        std::cout << "\n[FILTER] All " << num_queries << " queries have ground truth.\n";
    }
    time_truth = time(NULL);
    
    seqan::Index<seqan::Dna5String, seqan::FMIndex<>> fm_index;
    seqan::Dna5String ref_seq;
    if (method == "sae"){
        // ----- FMindex -----
        cout << "Building FM-Index:" << endl;
        seqan::assign(ref_seq, ref);
        // 使用 Dna5String 作为文本类型，FMIndex<> 作为索引类型
        fm_index = seqan::Index<seqan::Dna5String, seqan::FMIndex<>>(ref_seq);
        try {
            seqan::indexRequire(fm_index, seqan::FibreSA());
            std::cout << "[INFO] FM-index constructed successfully and Suffix Array loaded.\n";
        } catch (const seqan::Exception &e) {
            std::cerr << "[ERROR] Failed to construct or load Suffix Array for FM-index: " << e.what() << "\n";
            return 1; // 或者采取其他错误处理措施
        }
    }

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
    vector<size_t> all_node_accesses(num_queries);
    vector<vector<double>> all_radii_vecs(num_queries);

    #pragma omp parallel for default(none) \
        shared(queries, truth_positions, mtree, maxDist, fm_index, ref_seq, seed_len, anchor_index, anchors, use_all, start_idx, end_idx, last_random, anchor_radius, ref, method, dist_counts, all_node_accesses, all_radii_vecs, std::cout, num_queries) \
        reduction(+:sumTP, sumFP, sumFN, sum_avg_dist, sum_max_dist)
    for (size_t i = 0; i < num_queries; ++i) {
        const auto &q = queries[i];
        const auto &truth_pos = truth_positions[i];
        
        // === metrics ===
        size_t count = 0;  // 保存 dist_list.size()
        const auto &truth_str = posVecToStrVec(truth_positions[i]);
        vector<string> cand_str;
        
        // --- 临时变量，用于存储此线程的结果 ---
        size_t node_access_local = 0;
        std::vector<double> radii_vec_local;
        
        if (method == "mtree"){
            auto [pos, access, radii_vec] = retrieveCandidates_mtree(mtree, queries[i], maxDist);
            cand_str = posVecToStrVec(pos);
            
            // 缓存结果，而不是直接打印
            node_access_local = access;
            radii_vec_local = std::move(radii_vec);

        } else if (method == "sae"){
            cand_str = posVecToStrVec(retrieveCandidates_sae(fm_index, ref_seq, queries[i], maxDist, seed_len));
        } else{
            cand_str = posVecToStrVec(retrieveCandidates_anchor(queries[i], anchor_index, anchors[0].seq.size(), maxDist, use_all, start_idx, end_idx, last_random, anchor_radius = anchor_radius, &count));
        };

        dist_counts[i] = count;

        int tp = 0;
        int fp = 0;
        int fn = 0;

        // 调用统一函数进行计算
        // 传入 truth_str, cand_str, 容错参数, 以及要填充的 tp, fp, fn 变量
        metrics::calculate_position_metrics(truth_str, cand_str, maxDist, tp, fp, fn);
        // int tp = metrics::countTP(truth_str, cand_str);
        // int fp = metrics::countFP(truth_str, cand_str);
        // int fn = metrics::countFN(truth_str, cand_str);
        auto [avg_dist, max_dist] = metrics::evaluateDistances(queries[i], cand_str, ref);

        sumTP += tp;
        sumFP += fp;
        sumFN += fn;
        sum_avg_dist += avg_dist;
        sum_max_dist += max_dist;

        // 将缓存的结果安全地存入主向量
        //    由于 i 是唯一的，这不需要锁
        all_node_accesses[i] = node_access_local;
        all_radii_vecs[i] = std::move(radii_vec_local);
    }
    time_query = time(NULL);

    // 5. 现在，在循环之后，安全地打印所有缓存的输出
    cout << "\n--- Begin Buffered Output ---\n";
    for (size_t i = 0; i < num_queries; ++i) {
        if (method == "mtree") {
            std::cout << "M-Tree search node accesses: " << all_node_accesses[i] << std::endl;
            for(double r : all_radii_vecs[i]) {
                std::cout << "LeafNode Radius: " << r << std::endl;
            }
        }
    }
    cout << "--- End Buffered Output ---\n";
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

    // 内存访问
    globalLogger.report();

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

    cout << "Original sequence length: " << records[0].seq.size() <<endl;


    return 0;
}
