#ifndef ANCHOR_GEN_H
#define ANCHOR_GEN_H

#include <vector>
#include <string>
#include "fasta_utils_seqan.hpp"

// 从序列中生成随机 k-mer
std::string random_anchor_from_seq(const std::string& seq, size_t k);

// 从多个 FASTA 记录生成 anchors
std::vector<FastaRecord> generate_anchors_from_fasta(
    const std::vector<FastaRecord> &records,
    size_t k,
    size_t num_anchors);

// 从单一参考序列生成 anchors
std::vector<FastaRecord> generate_anchors_from_seq(
    const std::string &ref,
    const std::string &ref_id,
    size_t k,
    size_t num_anchors);

// anchor 生成的主入口（命令行接口）
int anchor_gen_main(int argc, char** argv);

#endif // ANCHOR_GEN_H