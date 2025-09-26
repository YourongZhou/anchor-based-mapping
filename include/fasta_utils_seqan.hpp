#ifndef FASTA_UTILS_SEQAN_HPP
#define FASTA_UTILS_SEQAN_HPP

#include <seqan/seq_io.h>
#include <string>
#include <vector>
#include <random>
#include <stdexcept>
#include <iostream>

struct FastaRecord {
    std::string id;
    std::string seq;
};

// 读取 fasta 文件，返回 vector<FastaRecord>
inline std::vector<FastaRecord> read_fasta_seqan(const std::string &path) {
    seqan::SeqFileIn seqFileIn;
    if (!seqan::open(seqFileIn, path.c_str())) {
        throw std::runtime_error("Cannot open FASTA file: " + path);
    }

    std::vector<FastaRecord> records;
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    try {
        seqan::readRecords(ids, seqs, seqFileIn);
    } catch (...) {
        throw std::runtime_error("Failed to read FASTA records from: " + path);
    }

    for (size_t i = 0; i < seqan::length(ids); ++i) {
        records.push_back(FastaRecord{
            std::string(seqan::begin(ids[i]), seqan::end(ids[i])),
            std::string(seqan::begin(seqs[i]), seqan::end(seqs[i]))
        });
    }

    return records;
}

// 从一条序列中随机选取 anchor（长度 k）
inline std::string random_anchor_from_seq(const std::string &seq, size_t k) {
    if (seq.size() < k) {
        throw std::runtime_error("Sequence shorter than anchor length k");
    }
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, seq.size() - k);
    size_t start = dist(rng);
    return seq.substr(start, k);
}

// 从 fasta 库生成所有 anchors（每条序列选一个 anchor）
inline std::vector<FastaRecord> generate_anchors_from_fasta(const std::vector<FastaRecord> &records, size_t k) {
    std::vector<FastaRecord> anchors;
    for (const auto &r : records) {
        try {
            std::string anchor_seq = random_anchor_from_seq(r.seq, k);
            anchors.push_back({r.id + "_anchor", anchor_seq});
        } catch (std::runtime_error &e) {
            std::cerr << "Warning: " << e.what() << " Skipping record " << r.id << std::endl;
        }
    }
    return anchors;
}

// write_fasta using SeqAn2 SeqFileOut
inline void write_fasta_seqan(const std::string &path, const std::vector<FastaRecord> &recs) {
    try {
        seqan::SeqFileOut seqFileOut(path.c_str());
        for (const auto &r : recs) {
            seqan::CharString id = r.id;
            seqan::CharString seq = r.seq;
            writeRecord(seqFileOut, id, seq);
        }
    } catch (seqan::Exception const & e) {
        std::cerr << "ERROR: writing FASTA " << path << " : " << e.what() << std::endl;
    }
}

#endif // FASTA_UTILS_SEQAN_HPP
