#include "tools.h"


// 从 reference 随机生成 num_queries 条长度为 query_len 的子串（模拟 reads）
std::vector<std::string> simulate_queries(const std::string &ref,
                                                 int num_queries,
                                                 size_t query_len) {
    if (ref.size() < query_len) {
        throw std::runtime_error("Reference shorter than query length!");
    }

    std::vector<std::string> queries;
    queries.reserve(num_queries);

    // static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, ref.size() - query_len);

    for (int i = 0; i < num_queries; ++i) {
        size_t start = dist(rng);
        queries.push_back(ref.substr(start, query_len));
    }
    return queries;
}

// 用 SeqAn2 找到 query 在 ref 中的所有 exact-match 位置
std::vector<int> find_all_occurrences(const std::string &query,
                                             const std::string &ref) {
    using namespace seqan;
    typedef String<char> TText;
    typedef String<char> TPattern;

    TText text = ref;
    TPattern pattern = query;

    Finder<TText> finder(text);
    Pattern<TPattern, Horspool> patt(pattern);

    std::vector<int> positions;
    while (find(finder, patt)) {
        positions.push_back(position(finder));
    }
    return positions;
}

// 找到所有 query 在某个距离以内的 neighbor
std::vector<int> find_all_occurrences_approx(
    const std::string &query,
    const std::string &ref,
    int max_dist
) {
    std::vector<int> positions;
    int qlen = query.size();
    int rlen = ref.size();

    if (qlen > rlen) return positions;

    for (int i = 0; i <= rlen - qlen; i++) {
        std::string window = ref.substr(i, qlen);
        int d = levenshtein(query, window);
        if (d <= max_dist) {
            positions.push_back(i);
        }
    }

    return positions;
}

// 小工具：把 int positions 转成 string id （用于复用现有 string-based metrics）
std::vector<std::string> posVecToStrVec(const std::vector<int> &pos) {
    std::vector<std::string> out;
    out.reserve(pos.size());
    for (int p : pos) out.push_back(std::to_string(p));
    return out;
}

// ======================== 构建函数 ========================
void build_mtree_from_ref(const std::string &ref, int k, MTree &tree) {
    std::unordered_set<std::string> inserted;  // 记录已插入的序列
    size_t count = 0;

    for (size_t i = 0; i + k <= ref.size(); ++i) {
        std::string subseq = ref.substr(i, k);
        if (inserted.find(subseq) != inserted.end())
            continue; // 已插入则跳过

        Substring a{subseq, static_cast<int>(i)};
        tree.add(a);
        inserted.insert(subseq);

        count++;
        if (count % 1000 == 0)
            std::cout << "Inserted " << count << " unique anchors into MTree." << std::endl;
    }

    std::cout << "MTree construction completed. Total unique anchors: " << count << std::endl;
}

// ======================== 查询函数 ========================
std::tuple<vector<int>, size_t, std::vector<double>> retrieveCandidates_mtree(
    MTree &mtree,
    const std::string &query,
    int maxDist)
{
    std::vector<int> results;

    Substring query_anchor{query, -1}; // pos = -1 表示查询
    
    // 1. 调用函数，'result_struct' 是包含 .matches 和 .nodeAccesses 的结构体
    auto result_struct = mtree.get_nearest_by_range(query_anchor, maxDist);

    // 2. 遍历结构体内部的 .matches 向量
    for (const auto &res : result_struct.matches) {
        // 'res' 现在直接是 result_item (包含 .data 和 .distance)
        results.push_back(res.data.pos);
    }

    // 3. (可选) 获取节点访问次数
    size_t accesses = result_struct.nodeAccesses;
    std::cout << "M-Tree search node accesses: " << accesses << std::endl;

    std::vector<double> radii = std::move(result_struct.leafNodeRadii);

    return {results, accesses, radii};
}


std::vector<int> retrieveCandidates_sae(
    seqan::Index<seqan::Dna5String, seqan::FMIndex<>> &fm_index, // 显式类型
    const seqan::Dna5String &ref_seq, // 显式类型
    const std::string &query,
    int maxDist,
    int seed_len) 
{
    cout << "Getting seed for " << query;
    // 显式使用 seqan::Score 和 seqan::Simple
    seqan::Score<int, seqan::Simple> scoring(1, -1, -1); // match=+1, mismatch=-1, gap=-1
    int xDropThreshold = maxDist * 4; // 经验值
    // ======================

    std::set<int> unique_positions;
    std::vector<int> results;

    std::string qseq_str = query;
    std::transform(qseq_str.begin(), qseq_str.end(), qseq_str.begin(), ::toupper);
    
    // 显式使用 seqan::Dna5String 和 seqan::assign
    seqan::Dna5String qseq;
    seqan::assign(qseq, qseq_str);

    // 显式类型定义
    typedef seqan::Seed<seqan::Simple> TSeed;
    
    // 显式使用 seqan::Finder
    seqan::Finder<seqan::Index<seqan::Dna5String, seqan::FMIndex<>>> finder(fm_index);

    // 遍历所有可能的种子

    // 计算种子数量（根据 maxDist，代表容错能力）
    //maxDist 越大，种子数量越多；容错能力越强
    seed_len = qseq_str.size() / maxDist;
    for (size_t i = 0; i + seed_len <= qseq_str.size(); i += seed_len) {
    // for (size_t i = 0; i + seed_len <= qseq_str.size(); i += 3) {
        globalLogger.accessIndex("seed");
        std::string seed_str = qseq_str.substr(i, seed_len);
        // cout << "seed: " << seed_str << " ";
        seqan::Dna5String seed;
        seqan::assign(seed, seed_str);
        
        seqan::clear(finder);
        
        // 显式使用 seqan::find, seqan::position
        while (seqan::find(finder, seed)) {
            cout << "found seed!";
            globalLogger.accessCandidate("extend");
            size_t seed_ini_pos_on_ref = seqan::position(finder);
            
            // 显式使用 TSeed 构造函数
            TSeed s(seed_ini_pos_on_ref, i, 
                    seed_ini_pos_on_ref + seed_len - 1, i + seed_len - 1);
            
            // 显式使用 seqan::extendSeed, seqan::EXTEND_BOTH, seqan::GappedXDrop
            seqan::extendSeed(s, ref_seq, qseq, seqan::EXTEND_BOTH, scoring, xDropThreshold, seqan::GappedXDrop());

            // 显式使用 seqan::beginPositionH
            unsigned sb = seqan::beginPositionH(s);
            
            // 存储起始位置
            unique_positions.insert((int)sb); 
        }
    }

    results.reserve(unique_positions.size());
    results.assign(unique_positions.begin(), unique_positions.end());
    
    return results;
}