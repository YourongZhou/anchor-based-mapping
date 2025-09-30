#include "fasta_utils_seqan.hpp"
#include "fastq_utils_seqan.hpp"
#include "levenshtein.hpp"
#include "metrics.h"

#include <iostream>
#include <vector>
#include <string>


// 计算 query -> anchors 距离
std::vector<std::pair<std::string, int>> queryAnchors(const std::string &query,
                                                      const std::vector<FastaRecord> &anchors) {
    std::vector<std::pair<std::string, int>> result;
    for (const auto &a : anchors) {
        int d = levenshtein(a.seq, query);
        result.push_back({a.id, d});
    }
    return result;
}

// 根据距离筛选 candidate
std::vector<std::string> retrieveCandidates(const std::vector<std::pair<std::string, int>> &anchorDist,
                                            int maxDist = 5) {
    std::vector<std::string> candidates;
    for (auto &p : anchorDist) {
        if (p.second <= maxDist) {
            candidates.push_back(p.first);
        }
    }
    return candidates;
}

int main() {
    // ===== Step 1: 读取 reference fasta =====
    std::string fasta_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fasta/ecoli.fa";
    size_t anchor_len = 150;
    int num_anchors = 10000;

    std::vector<FastaRecord> records;
    try {
        records = read_fasta_seqan(fasta_path);
    } catch (std::runtime_error &e) {
        std::cerr << "Error reading FASTA: " << e.what() << std::endl;
        return 1;
    }

    if (records.empty()) {
        std::cerr << "No records in " << fasta_path << std::endl;
        return 1;
    }

    // ===== Step 2: 生成 anchors =====
    auto anchors = generate_anchors_from_fasta(records, anchor_len, num_anchors);
    std::cout << "Generated " << anchors.size() << " anchors\n";
    // std::cout 

    // ===== Step 3: 读取 reads.fastq =====
    std::string fastq_path = "/home/luting/nfs/luting_data/AnchorBasedMapping/ecoli/fastq/ERR15404863_1.fastq";
    std::vector<FastqRecord> reads;
    try {
        reads = read_fastq_seqan(fastq_path);
    } catch (std::runtime_error &e) {
        std::cerr << "Error reading FASTQ: " << e.what() << std::endl;
        return 1;
    }

    if (reads.empty()) {
        std::cerr << "No reads in " << fastq_path << std::endl;
        return 1;
    }

    int totalTP = 0, totalFP = 0, totalFN = 0;

    // ===== Step 4: 对每条 read 做 mapping =====
    for (const auto &read : reads) {
        std::string query = read.seq.substr(0, anchor_len); // 截取前150bp
        auto anchorDist = queryAnchors(query, anchors);
        auto candidates = retrieveCandidates(anchorDist, 5);

        // 假设 ground truth 是第一个 anchor 对应的 read
        std::vector<std::string> truth = {anchors[0].id};

        // std::cout << "Read " << read.seq << " -> ";
        metrics::report(truth, candidates);

        totalTP += metrics::countTP(truth, candidates);
        totalFP += metrics::countFP(truth, candidates);
        totalFN += metrics::countFN(truth, candidates);
    }

    std::cout << "\n===== Overall Metrics =====\n";
    std::cout << "Total TP: " << totalTP << "\n";
    std::cout << "Total FP: " << totalFP << "\n";
    std::cout << "Total FN: " << totalFN << "\n";

    double recall = totalTP / double(totalTP + totalFN);
    double precision = totalTP / double(totalTP + totalFP);
    std::cout << "Overall Recall: " << recall << "\n";
    std::cout << "Overall Precision: " << precision << "\n";

    return 0;
}
