// src/part3.cpp
#include "part3.h"
#include <seqan/seq_io.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>

// 使用在 part2 中相同的 levenshtein
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

int retrieve_candidates(const char* anchor_map_path, const char* reference_fasta, const char* out_candidates_path) {
    try {
        // 读 anchor map
        std::ifstream fin(anchor_map_path);
        if (!fin) {
            std::cerr << "Cannot open anchor map: " << anchor_map_path << "\n";
            return 1;
        }
        std::vector<std::pair<std::string,int>> anchor_bounds;
        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string aid; int bound;
            if (!(iss >> aid >> bound)) continue;
            anchor_bounds.emplace_back(aid, bound);
        }
        fin.close();
        if (anchor_bounds.empty()) {
            std::cerr << "No anchors loaded from map\n";
            return 2;
        }

        // 读 reference fasta
        seqan::SeqFileIn seqFileIn;
        if (!seqan::open(seqFileIn, reference_fasta)) {
            std::cerr << "Cannot open reference fasta: " << reference_fasta << "\n";
            return 3;
        }

        std::ostream* out = &std::cout;
        std::ofstream fout;
        if (out_candidates_path && std::string(out_candidates_path).size()>0) {
            fout.open(out_candidates_path);
            if (!fout) {
                std::cerr << "Cannot open out candidates path " << out_candidates_path << "\n";
                return 4;
            }
            out = &fout;
        }

        // 为了减少重复计算：先把 anchors 的序列载入内存（anchor_id -> seq）
        // 这里我们需要 anchor 序列本身。简单做法：假设 anchors exist as a fasta file with same id names.
        // 为矬矬版：我们也期望在同目录下找到 anchors.fasta 与 anchor_map_path 同名的 anchors.fasta
        // 为简化：我们 attempt to open "anchors.fasta" in cwd. 如果需要更灵活可以把 anchor fasta 路径作为参数。
        std::unordered_map<std::string,std::string> anchor_seq_map;
        {
            const char* fallback_anchor_fasta = "anchors.fasta"; // 约定
            seqan::SeqFileIn anchorIn;
            if (!seqan::open(anchorIn, fallback_anchor_fasta)) {
                std::cerr << "Warning: cannot open fallback anchors.fasta to load anchor sequences. "
                          << "Distance checks will be approximate.\n";
                // We cannot compute distances without anchor sequences; return error for now.
                return 5;
            }
            while (!seqan::atEnd(anchorIn)) {
                seqan::CharString id; seqan::IupacString seq;
                seqan::readRecord(id, seq, anchorIn);
                anchor_seq_map[seqan::toCString(id)] = seqan::toCString(seq);
            }
        }

        // Iterate references
        while (!seqan::atEnd(seqFileIn)) {
            seqan::CharString rid; seqan::IupacString rseq;
            seqan::readRecord(rid, rseq, seqFileIn);
            std::string rid_s = seqan::toCString(rid);
            std::string rseq_s = seqan::toCString(rseq);
            bool inside = true;
            for (auto &p : anchor_bounds) {
                auto it = anchor_seq_map.find(p.first);
                if (it == anchor_seq_map.end()) {
                    std::cerr << "Anchor sequence " << p.first << " not found in anchors.fasta\n";
                    inside = false;
                    break;
                }
                int d = levenshtein(rseq_s, it->second);
                if (d > p.second) { inside = false; break; }
            }
            if (inside) {
                (*out) << rid_s << "\n";
            }
        }

        if (fout.is_open()) fout.close();
        return 0;
    } catch (std::exception &e) {
        std::cerr << "Exception in retrieve_candidates: " << e.what() << "\n";
        return 99;
    }
}
