#ifndef RETRIEVERS_H
#define RETRIEVERS_H

#include "tools.h"
#include "retriever_interface.h"
#include "query_grouper.h"
#include <seqan/index.h>
#include <vector>
#include <omp.h>

// M-Tree 检索器实现
class MTreeRetriever : public ICandidateRetriever {
public:
    MTreeRetriever(int minCap, int maxCap, double leafThresh, int compactMin, double compactFactor, int anchorLen, bool enableGrouping = false)
        : mtree(minCap, maxCap, leafThresh, compactMin, compactFactor, Distance(), SplitStrategyType()),
          anchor_len(anchorLen),
          enable_grouping(enableGrouping),
          routing_anchors_cached(false) {}

    void buildIndex(const std::string& reference) override {
        build_mtree_from_ref(reference, anchor_len, mtree);
        // 如果启用分组，预计算路由锚点
        if (enable_grouping) {
            routing_anchors = mtree.getRoutingAnchors();
            routing_anchors_cached = true;
            std::cout << "[MTreeRetriever] Routing anchors prepared: " << routing_anchors.size() << " anchors\n";
        }
    }

    CandidateResults retrieve(const std::string& query, int maxDist, bool requireDistance) override {
        auto [pos, access, radii_vec] = retrieveCandidates_mtree(mtree, query, maxDist);
        
        // 注意：这里不再设置固定的 0.0，而是让 Manager 通过 metrics::evaluateDistances 计算
        return {pos, 0.0, access, radii_vec};
    }

    /**
     * @brief 批量查询接口（Cache-Aware WORM 策略）
     * @param queries 查询字符串向量
     * @param maxDist 最大距离
     * @param requireDistance 是否需要计算距离
     * @return 查询结果向量，顺序与 queries 一致
     */
    std::vector<CandidateResults> retrieveBatch(
        const std::vector<std::string>& queries, 
        int maxDist, 
        bool requireDistance
    ) {
        std::vector<CandidateResults> results(queries.size());

        if (!enable_grouping || routing_anchors.empty()) {
            // 回退到原始并行模式
            #pragma omp parallel for schedule(dynamic)
            for (size_t i = 0; i < queries.size(); ++i) {
                results[i] = retrieve(queries[i], maxDist, requireDistance);
            }
            return results;
        }

        // Cache-Aware WORM 策略：先分桶，再按桶处理
        Distance distFunc;
        auto batches = QueryGrouper::groupQueries(queries, routing_anchors, distFunc);

        std::cout << "[MTreeRetriever] Grouped " << queries.size() << " queries into " 
                  << batches.size() << " buckets\n";

        // 按桶顺序处理（串行处理桶，桶内并行）
        for (const auto& batch : batches) {
            if (batch.query_indices.empty()) continue;

            // Hint: 这里未来可以加 prefetch 指令，把 anchors[batch.anchor_index] 对应的子树预取到缓存
            
            // 处理当前 Bucket 里的所有 Queries
            // 这些 Queries 极大概率会访问 M-Tree 的同一个分支 -> Cache Hit!
            #pragma omp parallel for schedule(static)
            for (size_t k = 0; k < batch.query_indices.size(); ++k) {
                int original_idx = batch.query_indices[k];
                const auto& q_str = queries[original_idx];
                
                // 调用原始的 M-Tree 查询逻辑
                // mtree 会自动利用 cache 里的热数据
                auto [pos, access, radii_vec] = retrieveCandidates_mtree(mtree, q_str, maxDist);
                results[original_idx] = {pos, 0.0, access, radii_vec};
            }
        }

        return results;
    }

    // reportStats 在并行环境下难以直接报告"最后一次"，建议由 Manager 统一处理聚合信息
    void reportStats() const override {}

    // 获取是否启用分组
    bool isGroupingEnabled() const { return enable_grouping; }

private:
    MTree mtree;
    int anchor_len;
    bool enable_grouping;
    std::vector<Data> routing_anchors;
    bool routing_anchors_cached;
};

// Anchor 检索器实现
class AnchorRetriever : public ICandidateRetriever {
public:
    AnchorRetriever(int numAnchors, int anchorLen, bool useAll, double startIdx, double endIdx, bool lastRandom, int anchorRadius)
        : num_anchors(numAnchors), anchor_len(anchorLen), use_all(useAll), 
          start_idx(startIdx), end_idx(endIdx), last_random(lastRandom), anchor_radius(anchorRadius) {}

    void buildIndex(const std::string& reference) override {
        ref_seq = reference;
        auto anchors_records = generate_anchors_from_seq(reference, "ref1", anchor_len, num_anchors);
        anchors = anchors_records;
        anchor_index = build_anchor_index(anchors_records, reference, anchor_len);
    }

    CandidateResults retrieve(const std::string& query, int maxDist, bool requireDistance) override {
        size_t dist_count = 0;
        seqan::Dna5String s_ref(ref_seq);
        return retrieveCandidates_anchor(query, anchor_index, s_ref, anchor_len, maxDist, 
                                        use_all, start_idx, end_idx, last_random, anchor_radius, &dist_count, requireDistance);
    }

private:
    int num_anchors;
    int anchor_len;
    bool use_all;
    double start_idx, end_idx;
    bool last_random;
    int anchor_radius;
    std::string ref_seq;
    std::vector<FastaRecord> anchors;
    std::unordered_map<std::string, std::vector<std::pair<int,int>>> anchor_index;
};

#include <mutex>

// SAE (FM-Index) 检索器实现
class SAERetriever : public ICandidateRetriever {
public:
    SAERetriever(int seedLen) : seed_len(seedLen) {}

    void buildIndex(const std::string& reference) override {
        seqan::assign(ref_seq, reference);
        fm_index = seqan::Index<seqan::Dna5String, seqan::FMIndex<>>(ref_seq);
        seqan::indexRequire(fm_index, seqan::FibreSA());
    }

    CandidateResults retrieve(const std::string& query, int maxDist, bool requireDistance) override {
        return retrieveCandidates_sae(fm_index, ref_seq, query, maxDist, seed_len, requireDistance);
    }

private:
    int seed_len;
    seqan::Dna5String ref_seq;
    seqan::Index<seqan::Dna5String, seqan::FMIndex<>> fm_index;
};

#endif // RETRIEVERS_H
