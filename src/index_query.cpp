// src/index_query.cpp
// Given anchors.fa and a query string (single read), compute anchor->distance map (JSON-like)
#include <iostream>
#include <map>
#include <string>
#include "fasta_utils_seqan.hpp"
#include "levenshtein.hpp"

// Usage: index_query <anchors.fa> <query_string>
// Output: JSON-like {"anchor_id":dist, ...}
int index_query_main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: index_query <anchors.fa> <query_string>\n";
        return 1;
    }
    std::string anchors_path = argv[1];
    std::string query = argv[2];

    auto anchors = read_fasta_seqan(anchors_path);
    if (anchors.empty()) {
        std::cerr << "No anchors read from " << anchors_path << std::endl;
        return 1;
    }

    bool first = true;
    std::cout << "{";
    for (const auto &a : anchors) {
        int d = levenshtein(a.seq, query);
        if (!first) std::cout << ", ";
        // minimal escape for quotes in id
        std::string id = a.id;
        for (char &c : id) if (c == '"') c = '\'';
        std::cout << "\"" << id << "\":" << d;
        first = false;
    }
    std::cout << "}\n";
    return 0;
}
