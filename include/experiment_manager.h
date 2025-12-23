#ifndef EXPERIMENT_MANAGER_H
#define EXPERIMENT_MANAGER_H

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "retriever_interface.h"
#include "retrievers.h"
#include "metrics.h"
#include "tools.h"
#include "access_log.hpp"

// 实验配置结构体
struct ExperimentConfig {
    std::string fasta_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fasta/ecoli.fa";
    std::string fastq_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fastq/ERR15404863_1.fastq";
    size_t truncate_ref_len = 10000;
    size_t anchor_len = 20;
    int num_anchors = 100;
    int num_queries = 100;
    std::vector<int> maxDist_list = {3};
    bool use_all = true;
    double start_idx = 0.0;
    double end_idx = 0.2;
    bool last_random = false;
    unsigned int seed = 42;
    int anchor_radius = 3;
    int min_node_capacity = 50;
    int max_node_capacity = -1;
    double leaf_radius_threshold = 0.0;
    int compactness_min_capacity = 25;
    double compactness_radius_factor = 1.1;
    std::string split_strategy_name = "mtree";
    bool auto_gen = false;
    int seed_len = 10;
    std::string method = "mtree";
    std::string index_dir = "./index_cache";
    bool force_rebuild = false;
    bool require_distance = false;
};

// 结果收集器
class MetricsCollector {
public:
    struct QueryResult {
        int tp, fp, fn;
        double avg_dist, max_dist;
        size_t node_access;
        std::vector<double> leaf_radii;
    };

    void addResult(const QueryResult& res) {
        #pragma omp critical
        {
            results.push_back(res);
            sumTP += res.tp;
            sumFP += res.fp;
            sumFN += res.fn;
            sumAvgDist += res.avg_dist;
            sumMaxDist += res.max_dist;
            sumNodeAccess += res.node_access;
        }
    }

    void report(int num_queries, const std::string& method) {
        if (num_queries <= 0) {
            std::cout << "\n===== Average per-query Metrics (" << method << ") =====\n";
            std::cout << "No queries to report.\n";
            return;
        }
        double recall = (sumTP + sumFN > 0) ? double(sumTP) / double(sumTP + sumFN) : 0.0;
        double precision = (sumTP + sumFP > 0) ? double(sumTP) / double(sumTP + sumFP) : 0.0;
        double fp_over_tp = (sumTP > 0) ? double(sumFP) / double(sumTP) : 0.0;

        std::cout << "\n===== Average per-query Metrics (" << method << ") =====\n";
        std::cout << "Average TP: " << (double)sumTP / num_queries << "\n";
        std::cout << "Average FP: " << (double)sumFP / num_queries << "\n";
        std::cout << "Average FN: " << (double)sumFN / num_queries << "\n";
        std::cout << "Average Recall: " << recall << "\n";
        std::cout << "Average Precision: " << precision << "\n";
        std::cout << "Average FP/TP: " << fp_over_tp << "\n";
        std::cout << "Average average distance: " << sumAvgDist / num_queries << "\n";
        std::cout << "Average maximum distance: " << sumMaxDist / num_queries << "\n";
        if (sumNodeAccess > 0) {
            std::cout << "Average M-Tree node access: " << (double)sumNodeAccess / num_queries << "\n";
        }

        std::cout << "\n===== All Individual Candidate Distances =====\n";
        std::cout << "DISTANCES_START:"; 
        for (const auto& res : results) {
            std::cout << " " << res.avg_dist;
        }
        std::cout << "\nDISTANCES_END\n";
    }

    void reset() {
        results.clear();
        sumTP = sumFP = sumFN = 0;
        sumAvgDist = sumMaxDist = 0.0;
        sumNodeAccess = 0;
    }

private:
    std::vector<QueryResult> results;
    int sumTP = 0, sumFP = 0, sumFN = 0;
    double sumAvgDist = 0.0, sumMaxDist = 0.0;
    size_t sumNodeAccess = 0;
};

// 实验管理器
class ExperimentManager {
public:
    ExperimentManager(const ExperimentConfig& config) : config(config) {
        retriever = createRetriever();
    }

    void printParams() {
        std::cout << "--- Paths ---\n";
        std::cout << "fasta_path: " << config.fasta_path << "\n";
        std::cout << "fastq_path: " << config.fastq_path << "\n";
        std::cout << "truncate_ref_len: " << config.truncate_ref_len << "\n";

        std::cout << "--- Anchor Settings ---\n";
        std::cout << "anchor_len: " << config.anchor_len << "\n";
        std::cout << "num_anchors: " << config.num_anchors << "\n";
        std::cout << "method: " << config.method << "\n";
    }

    void run() {
        // 0. 设置种子
        rng.seed(config.seed);
        
        // 打印参数确认 (可选，为了保持与原代码输出一致)
        printParams();

        // 1. 加载参考序列
        loadReference();

        // 2. 构建索引
        std::cout << "Building index using method: " << config.method << "..." << std::endl;
        retriever->buildIndex(reference);

        // 3. 运行批次查询
        for (int maxDist : config.maxDist_list) {
            runBatch(maxDist);
        }
    }

private:
    ExperimentConfig config;
    std::string reference;
    std::unique_ptr<ICandidateRetriever> retriever;
    MetricsCollector collector;

    std::unique_ptr<ICandidateRetriever> createRetriever() {
        if (config.method == "mtree") {
            return std::make_unique<MTreeRetriever>(
                config.min_node_capacity, config.max_node_capacity,
                config.leaf_radius_threshold, config.compactness_min_capacity,
                config.compactness_radius_factor, config.anchor_len);
        } else if (config.method == "anchor") {
            return std::make_unique<AnchorRetriever>(
                config.num_anchors, config.anchor_len, config.use_all,
                config.start_idx, config.end_idx, config.last_random, config.anchor_radius);
        } else if (config.method == "sae") {
            return std::make_unique<SAERetriever>(config.seed_len);
        }
        return nullptr;
    }

    void loadReference() {
        auto records = read_fasta_seqan(config.fasta_path);
        if (records.empty()) throw std::runtime_error("Empty FASTA file");
        reference = records[0].seq.substr(0, std::min(records[0].seq.size(), config.truncate_ref_len));
        std::cout << "Reference loaded, length: " << reference.size() << std::endl;
    }

    void runBatch(int maxDist) {
        std::cout << "---------------------" << "\n";
        std::cout << "--- START_BATCH_RUN ---" << "\n";
        std::cout << "---------------------" << "\n";
        std::cout << "query_maxDist: " << maxDist << "\n";

        collector.reset();
        globalLogger.reset(); // 确保每个 batch 的访问统计独立

        // 1. 生成/读取查询
        std::vector<std::string> queries;
        if (config.auto_gen) {
            queries = simulate_queries(reference, config.num_queries, config.anchor_len);
            std::cout << "Simulated " << queries.size() << " queries\n";
        } else {
            try {
                auto fastq_records = read_fastq_seqan(config.fastq_path);
                for (size_t i = 0; i < std::min((size_t)config.num_queries, fastq_records.size()); ++i) {
                    queries.push_back(fastq_records[fastq_records.size() - 1 - i].seq.substr(0, config.anchor_len));
                }
                std::cout << "Read " << queries.size() << " queries from FASTQ file: " << config.fastq_path << "\n";
            } catch (const std::exception& e) {
                std::cerr << "Error reading FASTQ: " << e.what() << "\n";
                return;
            }
        }

        // 2. 生成 Ground Truth
        std::cout << "Generating ground truths:\n";
        std::cout << "Number of queries: " << queries.size() << std::endl;
        std::vector<std::vector<int>> truth_positions(queries.size());
        
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < queries.size(); ++i) {
            truth_positions[i] = find_all_occurrences_approx(queries[i], reference, maxDist);
        }

        // 3. 过滤掉没有 Ground Truth 的查询 (保持与原逻辑一致)
        std::vector<std::string> filtered_queries;
        std::vector<std::vector<int>> filtered_truth;
        for (size_t i = 0; i < queries.size(); ++i) {
            if (!truth_positions[i].empty()) {
                filtered_queries.push_back(std::move(queries[i]));
                filtered_truth.push_back(std::move(truth_positions[i]));
            }
        }
        
        if (filtered_queries.size() < queries.size()) {
            std::cout << "\n[FILTER] Removed " << (queries.size() - filtered_queries.size()) 
                      << " queries with no ground truth (maxDist=" << maxDist << ").\n";
            queries = std::move(filtered_queries);
            truth_positions = std::move(filtered_truth);
        } else {
            std::cout << "\n[FILTER] All " << queries.size() << " queries have ground truth.\n";
        }

        // 4. 执行查询评估
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < queries.size(); ++i) {
            auto cand_results = retriever->retrieve(queries[i], maxDist, config.require_distance);
            
            const auto& truth_pos = truth_positions[i];
            const auto& cand_pos = cand_results.positions;
            
            int tp = 0, fp = 0, fn = 0;
            metrics::calculate_position_metrics(truth_pos, cand_pos, maxDist, tp, fp, fn);
            auto [avg_dist, max_dist] = metrics::evaluateDistances(queries[i], cand_pos, reference);

            collector.addResult({tp, fp, fn, avg_dist, (double)max_dist, cand_results.node_access, cand_results.leaf_node_radii});
        }

        // 5. 打印汇总报告
        collector.report(queries.size(), config.method);
        retriever->reportStats();
        globalLogger.report();
        std::cout << "---------------------" << "\n";
        std::cout << "--- END_BATCH_RUN ---" << "\n";
        std::cout << "---------------------" << "\n";
    }
};

#endif // EXPERIMENT_MANAGER_H

