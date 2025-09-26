#ifndef FASTQ_UTILS_SEQAN_HPP
#define FASTQ_UTILS_SEQAN_HPP

#include <seqan/seq_io.h>
#include <string>
#include <vector>
#include <stdexcept>

// 简单的 FastqRecord 结构体
struct FastqRecord {
    std::string id;
    std::string seq;
    std::string plus;
    std::string qual;
};

// 读取 FASTQ 文件，返回 vector<FastqRecord>
inline std::vector<FastqRecord> read_fastq_seqan(const std::string &path) {
    seqan::SeqFileIn seqFileIn;
    if (!seqan::open(seqFileIn, path.c_str())) {
        throw std::runtime_error("Cannot open FASTQ file: " + path);
    }

    std::vector<FastqRecord> records;
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString qual;

    try {
        while (!seqan::atEnd(seqFileIn)) {
            seqan::readRecord(id, seq, qual, seqFileIn);
            records.push_back(FastqRecord{
                std::string(seqan::begin(id), seqan::end(id)),
                std::string(seqan::begin(seq), seqan::end(seq)),
                "+", // plus 行在 FASTQ 格式中通常是一个 "+"
                std::string(seqan::begin(qual), seqan::end(qual))
            });
        }
    } catch (...) {
        throw std::runtime_error("Failed to read FASTQ records from: " + path);
    }

    return records;
}

#endif // FASTQ_UTILS_SEQAN_HPP
