// src/part2.cpp
#include "part2.h"
#include <seqan/seq_io.h>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

// 简单 Levenshtein edit distance, 两个 std::string
static int levenshtein(const std::string &a, const std::string &b) {
    size_t n = a.size(), m = b.size();
    if (n == 0) return (int)m;
    if (m == 0) return (int)n;
    std::vector<int> prev(m+1), cur(m+1);
    for (size_t j=0;j<=m;++j) prev[j] = j;
    for (size_t i=1;i<=n;++i) {
        cur[0] = i;
        for (size_t j=1;j<=m;++j) {
            int cost = (a[i-1]==b[j-1])?0:1;
            cur[j] = std::min({ prev[j] + 1, cur[j-1] + 1, prev[j-1] + cost });
        }
        std::swap(prev, cur);
    }
    return prev[m];
}

int compute_anchor_distances(const char* query_seq, const char* anchor_fasta, const char* out_map_path) {
    try {
        // load anchors
        seqan::SeqFileIn seqFileIn;
        if (!seqan::open(seqFileIn, anchor_fasta)) {
            std::cerr << "Cannot open anchor fasta: " << anchor_fasta << "\n";
            return 1;
        }

        std::vector<std::pair<std::string,std::string>> anchors;
        while (!seqan::atEnd(seqFileIn)) {
            seqan::CharString id; seqan::IupacString seq;
            seqan::readRecord(id, seq, seqFileIn);
            anchors.emplace_back(std::string(seqan::toCString(id)), std::string(seqan::toCString(seq)));
        }

        std::ostream* out = &std::cout;
        std::ofstream fout;
        if (out_map_path && std::string(out_map_path).size()>0) {
            fout.open(out_map_path);
            if (!fout) {
                std::cerr << "Cannot open out map path " << out_map_path << "\n";
                return 2;
            }
            out = &fout;
        }

        std::string q(query_seq);
        for (auto &p : anchors) {
            int d = levenshtein(q, p.second);
            (*out) << p.first << "\t" << d << "\n";
        }
        if (fout.is_open()) fout.close();
        return 0;
    } catch (std::exception &e) {
        std::cerr << "Exception in compute_anchor_distances: " << e.what() << "\n";
        return 99;
    }
}
