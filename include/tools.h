#ifndef TOOLS_H
#define TOOLS_H

#include "rng.h"
#include "fasta_utils_seqan.hpp"
#include "fastq_utils_seqan.hpp"
#include "levenshtein.hpp"
#include "metrics.h"
#include "anchor_gen.h"
#include "candidate_retriever.h"
#include "index_query.h"
#include "mtree.h"
#include "functions.h"
#include "access_log.hpp"
#include "data_types.h"

// #define SEQAN_NO_INCLUDE_OMP
#include <seqan/find.h> 
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/align.h>
#include <seqan/store.h>
#include <seqan/stream.h>

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <nlohmann/json.hpp>

// 从 reference 随机生成 num_queries 条长度为 query_len 的子串（模拟 reads）
std::vector<std::string> simulate_queries(const std::string &ref,
                                                 int num_queries,
                                                 size_t query_len);

// 用 SeqAn2 找到 query 在 ref 中的所有 exact-match 位置
std::vector<int> find_all_occurrences(const std::string &query,
                                             const std::string &ref);

// 找到所有 query 在某个距离以内的 neighbor
std::vector<int> find_all_occurrences_approx(
    const std::string &query,
    const std::string &ref,
    int max_dist
);

// 小工具：把 int positions 转成 string id （用于复用现有 string-based metrics）
std::vector<std::string> posVecToStrVec(const std::vector<int> &pos);

// ======================== 构建函数 ========================
void build_mtree_from_ref(const std::string &ref, int k, MTree &tree);

// ======================== 查询函数 ========================
std::tuple<std::vector<int>, size_t, std::vector<double>> retrieveCandidates_mtree(
    MTree &mtree,
    const std::string &query,
    int maxDist);


CandidateResults retrieveCandidates_sae(
    seqan::Index<seqan::Dna5String, seqan::FMIndex<>> &fm_index, // 显式类型
    const seqan::Dna5String &ref_seq, // 显式类型
    const std::string &query,
    int maxDist,
    int seed_len,
    bool require_distance);

std::string generate_index_filename(
    const std::string& method,
    const std::string& index_dir,
    size_t ref_len,
    size_t anchor_len,
    int num_anchors);

// 辅助函数，检查文件是否存在
#include <sys/stat.h>
bool file_exists(const std::string& filename);

std::string get_fm_index_cache_path(size_t ref_len);

// 写入缓存的函数
bool write_fm_index_cache(const seqan::Index<seqan::Dna5String, seqan::FMIndex<>>& fm_index, size_t ref_len);

// 读取缓存的函数
bool load_fm_index_cache(seqan::Index<seqan::Dna5String, seqan::FMIndex<>>& fm_index, size_t ref_len);

// 辅助函数：将逗号分隔的字符串解析为整数列表
std::vector<int> parseIntList(const std::string& str);

#endif
