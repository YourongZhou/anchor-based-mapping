// src/candidate_retriever.cpp
// Input: anchors.fa, anchor_dist_json_file, refs.fasta
// Output: list of candidate reference IDs (stdout), one per line
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include "fasta_utils_seqan.hpp"
#include "levenshtein.hpp"

// crude JSON map parser expecting {"id":num,...}
static std::map<std::string,int> parse_simple_json_map(const std::string &s) {
    std::map<std::string,int> out;
    size_t i = 0, n = s.size();
    while (i < n && s[i] != '{') ++i;
    if (i == n) return out;
    ++i;
    while (i < n) {
        while (i<n && isspace((unsigned char)s[i])) ++i;
        if (i<n && s[i]=='}') break;
        if (i<n && s[i]=='"') {
            ++i;
            std::string id;
            while (i<n && s[i]!='"') id.push_back(s[i++]);
            ++i;
            while (i<n && (s[i]==':' || isspace((unsigned char)s[i]))) ++i;
            std::string num;
            while (i<n && (s[i]=='-' || isdigit((unsigned char)s[i]))) num.push_back(s[i++]);
            if (!num.empty()) {
                int val = std::stoi(num);
                out[id] = val;
            }
            while (i<n && s[i] != ',' && s[i] != '}') ++i;
            if (i<n && s[i]==',') ++i;
        } else ++i;
    }
    return out;
}

int candidate_retriever_main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: candidate_retriever <anchors.fa> <anchor_dist_json_file> <refs.fasta>\n";
        return 1;
    }
    std::string anchors_path = argv[1];
    std::string json_file = argv[2];
    std::string refs = argv[3];

    auto anchors = read_fasta_seqan(anchors_path);
    std::map<std::string,std::string> anchor_seq;
    for (const auto &a : anchors) anchor_seq[a.id] = a.seq;

    // read JSON file
    std::ifstream ifs(json_file);
    if (!ifs) {
        std::cerr << "Cannot open " << json_file << std::endl;
        return 1;
    }
    std::string js; std::ostringstream oss;
    oss << ifs.rdbuf();
    js = oss.str();
    auto amap = parse_simple_json_map(js);
    if (amap.empty()) {
        std::cerr << "No anchors parsed from " << json_file << std::endl;
        return 1;
    }

    // build vector of (anchor_seq, threshold)
    std::vector<std::pair<std::string,int>> queries;
    for (auto &kv : amap) {
        auto it = anchor_seq.find(kv.first);
        if (it == anchor_seq.end()) {
            std::cerr << "Warning: anchor id " << kv.first << " not found in anchors.fa\n";
            continue;
        }
        queries.emplace_back(it->second, kv.second);
    }
    if (queries.empty()) {
        std::cerr << "No valid anchors for query\n";
        return 1;
    }

    auto ref_records = read_fasta_seqan(refs);
    size_t cand_count = 0;
    for (const auto &r : ref_records) {
        bool ok = true;
        for (const auto &q : queries) {
            int d = min_edit_distance_window(q.first, r.seq);
            if (d > q.second) { ok = false; break; }
        }
        if (ok) {
            std::cout << r.id << "\n";
            ++cand_count;
        }
    }
    std::cerr << "Candidates found: " << cand_count << std::endl;
    return 0;
}
