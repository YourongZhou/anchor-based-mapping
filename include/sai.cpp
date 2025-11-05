/* 
===================================================================
Task: realise seed-and-index process

- Input: .fasta and .fastq, which represents reference and query respectively
- Process: 
    - FM-index
    - For each query:
        - seeding
        - match seeds
        - seed extension
        - record outcome (?)
- Analysis: comparison, ...

===================================================================
*/ 

#include "fasta_utils_seqan.hpp"
#include "fastq_utils_seqan.hpp"
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <iostream>
#include <algorithm>

using namespace seqan;

int main() {
    std::string fasta_path = "/home/minghao/project/Projects/Mapping/ecoli/fasta/ecoli.fa";   // reference 文件路径
    std::string fastq_path = "/home/minghao/project/Projects/Mapping/ecoli/fastq/ERR15404863_1.fastq";       // query 文件路径
    size_t truncate_ref_len = 100000;

    // Reading fasta
    std::vector<FastaRecord> ref_records;
    try {
        ref_records = read_fasta_seqan(fasta_path);
    } catch (std::runtime_error &e) {
        std::cerr << "Error reading FASTA: " << e.what() << "\n";
        return 1;
    }
    if (ref_records.empty()) {
        std::cerr << "No sequences found in FASTA file.\n";
        return 1;
    }

    // Making reference
    std::string ref = ref_records[0].seq.substr(0, std::min(ref_records[0].seq.size(), truncate_ref_len));
                       //std::min(ref_records[0].seq.size(), truncate_ref_len));
    std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);

    std::cout << "[INFO] Using reference sequence: " << ref_records[0].id << "\n";
    std::cout << "[INFO] Reference length (possibly truncated) = " << ref.size() << "\n";

    // Reading fastq
    std::vector<FastqRecord> query_records;
    try {
        query_records = read_fastq_seqan(fastq_path);
    } catch (std::runtime_error &e) {
        std::cerr << "Error reading FASTQ: " << e.what() << "\n";
        return 1;
    }
    if (query_records.empty()) {
        std::cerr << "No reads found in FASTQ file.\n";
        return 1;
    }

    std::cout << "[INFO] Loaded " << query_records.size() << " query reads.\n";
    std::cout << "[INFO] First query ID: " << query_records[0].id << "\n";
    std::cout << "[INFO] First query sequence (first 50 bp): "
              << query_records[0].seq.substr(0, std::min<size_t>(50, query_records[0].seq.size()))
              << "\n";

    //============================================================
    // Build FM-index for reference
    Dna5String ref_seq;
    assign(ref_seq, ref);
    Index<Dna5String, FMIndex<>> fm_index(ref_seq);

    indexRequire(fm_index, FibreSA());

    std::cout << "[INFO] FM-index constructed successfully. \n";

    //============================================================
    // Start matching
    size_t seed_len;
    size_t max_queries = 150;
    size_t processed_queries = 0;
    double match_threshold = 0.92; // Set threshold for successful match

    for(seed_len = 15; seed_len <= 15; seed_len++) { // In memory of the NBA legend Kobe Bryant
        size_t matched_count = 0;
        processed_queries = 0;

        for(const auto &q : query_records) {
            // Fetch a query and ++processed_queries. When this number reaches max, go to next seed length
            if(processed_queries >= max_queries) break;
            processed_queries++;

            typedef Seed<Simple> TSeed;
            SeedSet<TSeed> seedSet;

            // Slice the query
            for(size_t i = 0; i + seed_len <= q.seq.size(); i += 3) {
                std::string seed_str = q.seq.substr(i, seed_len);
                Dna5String seed;
                assign(seed, seed_str);

                // Find seed in FM-index
                Finder<Index<Dna5String, FMIndex<>>> finder(fm_index);
                clear(finder);

                while(find(finder, seed)) {
                    size_t seed_ini_pos_on_ref = position(finder);

                    // Create Seed object and add to seedSet
                    TSeed s(seed_ini_pos_on_ref, i, seed_ini_pos_on_ref + seed_len - 1,
                            i + seed_len - 1);

                            // Extend seed
                    Score<int, Simple> scoring(1, -1, -1);
                    extendSeed(s, ref_seq, q.seq, EXTEND_BOTH, scoring, 2, GappedXDrop());
                
                    if (endPositionH(s) > length(ref_seq) || endPositionV(s) > length(q.seq)) {
                        std::cerr << "[Warning] Seed out of range! Skipping..." << std::endl;
                        continue;
                    }
                    addSeed(seedSet, s, Single());
                    std::cout << "[INFO] SeedSet size: " << length(seedSet) << std::endl;
                }
            }
            
            // 1. 定义 SeedSet 的常量迭代器类型
                // 根据您的要求，明确指定 SeqAn 标准常量迭代器
                using TSeedSet = SeedSet<TSeed>;
                typedef typename seqan::Iterator<TSeedSet, seqan::Standard>::Type TIterator;

                // 如果 seedSet 是空的，则跳过 Chaining
                if (length(seedSet) == 0) {
                    std::cout << "--- SeedSet for query " << processed_queries << " is empty. Skipping chaining. ---" << std::endl;
                    continue;
                }
                // 2. 输出 SeedSet 的信息
                
                std::cout << "--- SeedSet Contents (Total: " << length(seedSet) << " Seeds) ---" << std::endl;

                TIterator it_end = seqan::end(seedSet, seqan::Standard());

                // 3. 显式 for 循环遍历
                for (TIterator it = seqan::begin(seedSet, seqan::Standard()); 
                    it != it_end; 
                    ++it)
                {
                    // 强制使用 const 引用来获取当前的种子（使用变量名 s）
                    // 这解决了 "discards qualifiers" 的编译错误
                    const TSeed& s = *it; 
                    
                    // 4. 输出种子的位置信息
                    std::cout << "Seed | [Ref: "
                            << seqan::beginPositionH(s) << " - " << seqan::endPositionH(s) // H = Horizontal (参考序列)
                            << "] | [Query: "
                            << seqan::beginPositionV(s) << " - " << seqan::endPositionV(s) // V = Vertical (查询序列)
                            << std::endl;
                }

            std::cout << "----------------------------------------------------" << std::endl;

            // Global chaining
            //seqan::simplify(seedSet);
            String<TSeed> chainResult;
            chainSeedsGlobally(chainResult, seedSet, SparseChaining());

            // Calculate coverage
            /*size_t chain_len = 0;
            for (auto &s : chainResult)
                chain_len += length(s);

            double coverage_ratio = static_cast<double>(chain_len) / q.seq.size();

            if (coverage_ratio >= match_threshold)
                matched_count++;*/
        }

        /*std::cout << "[INFO] Seed length: " << seed_len
                  << ", matched " << matched_count
                  << "/" << processed_queries
                  << " queries, recall: "
                  << static_cast<double>(matched_count) / processed_queries
                  << std::endl;*/
    }

    return 0;
}

//g++ seedandindex.cpp -o seedandindex
