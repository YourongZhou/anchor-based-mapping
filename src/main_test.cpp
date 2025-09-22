// src/main_test.cpp
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

// 声明外部库函数（动态链接或静态链接都可）
extern "C" {
    int generate_anchors(const char* input_fasta, const char* out_anchor_fasta, int num_anchors, unsigned int seed);
    int compute_anchor_distances(const char* query_seq, const char* anchor_fasta, const char* out_map_path);
    int retrieve_candidates(const char* anchor_map_path, const char* reference_fasta, const char* out_candidates_path);
}

std::map<std::string,int> read_anchor_map(const std::string &path) {
    std::map<std::string,int> m;
    std::ifstream fin(path);
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string id; int d;
        iss >> id >> d;
        m[id]=d;
    }
    return m;
}

std::set<std::string> read_candidates(const std::string &path) {
    std::set<std::string> s;
    std::ifstream fin(path);
    std::string line;
    while (std::getline(fin,line)) {
        if (line.empty()) continue;
        s.insert(line);
    }
    return s;
}

// 读取 truth mapping 文件，格式: query_id \t true_ref_id
std::map<std::string,std::string> read_truth(const std::string &path) {
    std::map<std::string,std::string> m;
    std::ifstream fin(path);
    std::string line;
    while (std::getline(fin,line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string q,t;
        iss >> q >> t;
        m[q]=t;
    }
    return m;
}

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <reference_fasta> <query_fasta> <num_anchors> <truth_map(optional:path or '-')> <seed>\n";
        return 1;
    }
    std::string ref_fasta = argv[1];
    std::string query_fasta = argv[2];
    int num_anchors = std::stoi(argv[3]);
    std::string truth_path = argv[4];
    unsigned int seed = std::stoul(argv[5]);

    // 1. generate anchors
    std::string anchors_fasta = "anchors.fasta";
    if (generate_anchors(ref_fasta.c_str(), anchors_fasta.c_str(), num_anchors, seed) != 0) {
        std::cerr << "generate_anchors failed\n"; return 2;
    }
    std::cout << "Anchors generated: " << anchors_fasta << "\n";

    // For each query in query_fasta: compute anchor distances -> write map file, then retrieve candidates.
    // We'll read queries using simple seqan reading here by calling part2 for each sequence.
    // For simplicity, parse query fasta quickly using seqan functions inline (so link seqan in this binary)
    // To keep this file short, we call compute_anchor_distances via cli output file, then retrieve candidates.

    // We'll iterate queries by reading the query fasta with seqan (so include seqan here)
    // But to avoid adding seqan headers here (already linked), we can just call external small program -
    // For clarity, use a very small parsing using shell: (but avoid shell). Simpler: re-use compute_anchor_distances
    // which expects sequence string, so we need to read sequences ourselves. Let's quickly read query fasta using seqan.
    #include <seqan/seq_io.h> // placed inside file for brevity (allowed)

    seqan::SeqFileIn qIn;
    if (!seqan::open(qIn, query_fasta.c_str())) {
        std::cerr << "Cannot open query fasta: " << query_fasta << "\n"; return 3;
    }

    size_t total_queries = 0;
    size_t TP = 0;
    size_t total_candidates = 0;
    std::map<std::string,std::string> truth;
    if (truth_path != "-" && truth_path.size()>0) {
        truth = read_truth(truth_path);
    }

    while (!seqan::atEnd(qIn)) {
        seqan::CharString qid; seqan::IupacString qseq;
        seqan::readRecord(qid, qseq, qIn);
        std::string qid_s = seqan::toCString(qid);
        std::string qseq_s = seqan::toCString(qseq);

        std::string map_out = "anchor_map.tmp.txt";
        if (compute_anchor_distances(qseq_s.c_str(), anchors_fasta.c_str(), map_out.c_str()) != 0) {
            std::cerr << "compute_anchor_distances failed for query " << qid_s << "\n"; continue;
        }

        std::string cand_out = "candidates.tmp.txt";
        if (retrieve_candidates(map_out.c_str(), ref_fasta.c_str(), cand_out.c_str()) != 0) {
            std::cerr << "retrieve_candidates failed for query " << qid_s << "\n"; continue;
        }

        auto candidates = read_candidates(cand_out);
        total_candidates += candidates.size();
        total_queries += 1;

        if (truth.count(qid_s)) {
            const std::string &true_ref = truth[qid_s];
            if (candidates.count(true_ref)) TP++;
        }

        // clean tmp files if you want
        remove(map_out.c_str());
        remove(cand_out.c_str());
    }

    std::cout << "Queries processed: " << total_queries << "\n";
    if (!truth.empty()) {
        std::cout << "TP (true ref found in candidates): " << TP << " / " << total_queries << "\n";
    } else {
        std::cout << "No truth provided; TP not computed.\n";
    }
    std::cout << "Average candidate size per query: " << (total_queries? (double)total_candidates/total_queries : 0.0) << "\n";

    return 0;
}
