// src/part1.cpp
#include "part1.h"
#include <seqan/seq_io.h>
#include <vector>
#include <random>
#include <iostream>

int generate_anchors(const char* input_fasta, const char* out_anchor_fasta, int num_anchors, unsigned int seed) {
    try {
        seqan::SeqFileIn seqFileIn;
        if (!seqan::open(seqFileIn, input_fasta)) {
            std::cerr << "Cannot open input FASTA: " << input_fasta << "\n";
            return 1;
        }

        std::vector<std::string> ids;
        std::vector<std::string> seqs;
        seqan::StringSet<seqan::CharString> seqIds;
        seqan::StringSet<seqan::IupacString> sequences;

        // read all records
        while (!seqan::atEnd(seqFileIn)) {
            seqan::CharString id; seqan::IupacString seq;
            seqan::readRecord(id, seq, seqFileIn);
            ids.push_back(seqan::toCString(id));
            seqs.push_back(seqan::toCString(seq));
        }

        if ((int)seqs.size() == 0) {
            std::cerr << "No sequences found in input FASTA\n";
            return 2;
        }

        std::mt19937 rng(seed);
        std::uniform_int_distribution<size_t> dist(0, seqs.size()-1);
        std::vector<size_t> chosen;
        chosen.reserve(num_anchors);

        // simple reservoir-like sampling without guarantee of uniqueness -> ensure unique picks
        std::unordered_set<size_t> chosen_set;
        while ((int)chosen_set.size() < num_anchors && chosen_set.size() < ids.size()) {
            chosen_set.insert(dist(rng));
        }

        // write anchors to out_anchor_fasta
        seqan::SeqFileOut seqFileOut;
        if (!seqan::open(seqFileOut, out_anchor_fasta)) {
            std::cerr << "Cannot open output anchor FASTA: " << out_anchor_fasta << "\n";
            return 3;
        }

        for (auto idx : chosen_set) {
            seqan::CharString id = seqan::CharString(ids[idx].c_str());
            seqan::IupacString seq = seqan::IupacString(seqs[idx].c_str());
            seqan::writeRecord(seqFileOut, id, seq);
        }

        return 0;
    } catch (std::exception &e) {
        std::cerr << "Exception in generate_anchors: " << e.what() << "\n";
        return 99;
    }
}
