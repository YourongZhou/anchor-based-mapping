#ifndef RETRIEVERS_H
#define RETRIEVERS_H

#include "tools.h"
#include "retriever_interface.h"
#include <seqan/index.h>

// M-Tree 检索器实现
class MTreeRetriever : public ICandidateRetriever {
public:
    MTreeRetriever(int minCap, int maxCap, double leafThresh, int compactMin, double compactFactor, int anchorLen)
        : mtree(minCap, maxCap, leafThresh, compactMin, compactFactor, Distance(), SplitStrategyType()),
          anchor_len(anchorLen) {}

    void buildIndex(const std::string& reference) override {
        build_mtree_from_ref(reference, anchor_len, mtree); 
    }

    CandidateResults retrieve(const std::string& query, int maxDist, bool requireDistance) override {
        auto [pos, access, radii_vec] = retrieveCandidates_mtree(mtree, query, maxDist);
        
        double avg_dist = 0;
        // 注意：这里不再修改成员变量，保持线程安全
        return {pos, avg_dist, access, radii_vec};
    }

    // reportStats 在并行环境下难以直接报告“最后一次”，建议由 Manager 统一处理聚合信息
    void reportStats() const override {}

private:
    MTree mtree;
    int anchor_len;
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
