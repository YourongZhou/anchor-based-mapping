// src/anchor_gen.cpp
// Generate anchors (random k-length substrings) from reference FASTA (SeqAn2 used for IO)
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <cstdlib>
#include <unordered_set>
#include "rng.h"
#include "fasta_utils_seqan.hpp"
#include "anchor_gen.h"

// 从 fasta 库生成所有 anchors
std::vector<FastaRecord> generate_anchors_from_fasta(
    const std::vector<FastaRecord> &records,
    size_t k,
    size_t num_anchors) 
{
    if (records.empty()) {
        throw std::runtime_error("No records provided for anchor generation");
    }

    std::vector<FastaRecord> anchors;
    anchors.reserve(num_anchors);

    std::unordered_set<std::string> seen;  // 用来去重
    seen.reserve(num_anchors * 2);         // 预留容量，避免频繁扩展

    std::uniform_int_distribution<size_t> record_dist(0, records.size() - 1);

    size_t generated = 0;
    size_t attempts = 0; // 防止死循环

    while (generated < num_anchors) {
        const auto &rec = records[record_dist(rng)];
        if (rec.seq.size() < k) {
            // 跳过太短的序列
            continue;
        }

        std::string anchor_seq = random_anchor_from_seq(rec.seq, k);

        // 如果重复，就跳过
        if (seen.find(anchor_seq) != seen.end()) {
            attempts++;
            if (attempts > num_anchors * 10) {
                throw std::runtime_error("Too many duplicate anchors, cannot reach desired num_anchors");
            }
            continue;
        }

        // 记录 & 保存
        seen.insert(anchor_seq);
        anchors.push_back({rec.id + "_anchor_" + std::to_string(generated), anchor_seq});
        generated++;
    }

    return anchors;
}


// 从单条 reference 序列直接生成 anchors
std::vector<FastaRecord> generate_anchors_from_seq(
    const std::string &ref,
    const std::string &ref_id,
    size_t k,
    size_t num_anchors)
{
    // 包装成一个临时的 FastaRecord
    std::vector<FastaRecord> tmp = {{ref_id, ref}};

    // 复用 generate_anchors_from_fasta
    return generate_anchors_from_fasta(tmp, k, num_anchors);
}


int anchor_gen_main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: anchor_gen <refs.fasta> <out_anchors.fa> <k> <anchors_per_seq>\n";
        return 1;
    }
    std::string refs = argv[1];
    std::string out = argv[2];
    int k = std::atoi(argv[3]);
    int anchors_per_seq = std::atoi(argv[4]);
    if (k <= 0 || anchors_per_seq <= 0) {
        std::cerr << "k and anchors_per_seq must be positive\n";
        return 1;
    }

    auto recs = read_fasta_seqan(refs);
    if (recs.empty()) {
        std::cerr << "No reference sequences read from " << refs << std::endl;
        return 1;
    }

    // std::mt19937_64 rng(123456789ULL);
    std::vector<FastaRecord> anchors;
    for (const auto &r : recs) {
        const std::string &s = r.seq;
        if ((int)s.size() < k) continue;
        std::uniform_int_distribution<int> dist(0, (int)s.size() - k);
        for (int i = 0; i < anchors_per_seq; ++i) {
            int off = dist(rng);
            std::string sub = s.substr(off, k);
            std::string id = r.id + "|off=" + std::to_string(off);
            anchors.push_back({id, sub});
        }
    }

    write_fasta_seqan(out, anchors);
    std::cerr << "Wrote " << anchors.size() << " anchors to " << out << std::endl;
    return 0;
}


// 从一条序列中随机选取 anchor（长度 k）
inline std::string random_anchor_from_seq(const std::string &seq, size_t k) {
    if (seq.size() < k) {
        throw std::runtime_error("Sequence shorter than anchor length k");
    }
    // static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, seq.size() - k);
    size_t start = dist(rng);
    return seq.substr(start, k);
}