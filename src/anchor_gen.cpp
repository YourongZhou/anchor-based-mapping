// src/anchor_gen.cpp
// Generate anchors (random k-length substrings) from reference FASTA (SeqAn2 used for IO)
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <cstdlib>
#include "fasta_utils_seqan.hpp"

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

    std::mt19937_64 rng(123456789ULL);
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
