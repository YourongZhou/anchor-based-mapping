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
    // --- 1. (并行) 阶段 1: 提取所有 k-mer 及其位置 ---
    // (注意：这里包含重复项，我们稍后处理)
    // 预先分配内存以避免并行时的 re-allocation
    std::vector<std::pair<std::string, int>> all_kmers;
    all_kmers.resize(ref.size() - k + 1);

    // 使用 OpenMP 并行化 substr 提取
    // 这是可以安全并行的，因为每个线程写入 all_kmers[i] 的不同位置
    #pragma omp parallel for default(none) shared(ref, k, all_kmers)
    for (size_t i = 0; i <= ref.size() - k; ++i) {
        all_kmers[i] = {ref.substr(i, k), static_cast<int>(i)};
    }

    // --- 2. (串行) 阶段 2: 消除重复并构建 M-Tree ---
    // 此时，我们回到单线程，因为 tree.add() 必须串行
    std::cout << "Parallel k-mer extraction complete. Deduplicating and building tree..." << std::endl;

    std::unordered_set<std::string> inserted; // 记录已插入的序列
    size_t count = 0;

    // 遍历我们收集的 k-mer (包含重复)
    for (const auto& kv_pair : all_kmers) {
        const std::string& subseq = kv_pair.first;
        int pos = kv_pair.second;

        // 检查是否已经插入过
        if (inserted.find(subseq) != inserted.end())
            continue; // 已插入则跳过

        // 这是该 k-mer 第一次出现，插入 M-Tree
        Substring a{subseq, pos};
        tree.add(a);
        inserted.insert(subseq); // 标记为已插入

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
        #pragma omp critical (GlobalLoggerLock)
        {
            globalLogger.accessIndex("seed");
        }
        std::string seed_str = qseq_str.substr(i, seed_len);
        // cout << "seed: " << seed_str << " ";
        seqan::Dna5String seed;
        seqan::assign(seed, seed_str);
        
        seqan::clear(finder);
        
        // 显式使用 seqan::find, seqan::position
        while (seqan::find(finder, seed)) {
            // cout << "found seed!";
            #pragma omp critical (GlobalLoggerLock)
            {
                globalLogger.accessCandidate("extend");
            }
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

std::string generate_index_filename(
    const std::string& method,
    const std::string& index_dir,
    size_t ref_len,
    size_t anchor_len,
    int num_anchors)
{
    std::stringstream ss;
    ss << index_dir << "/";
    ss << method << "_ref" << ref_len;
    
    if (method == "anchor") {
        ss << "_k" << anchor_len << "_n" << num_anchors;
    }
    
    ss << ".index";
    return ss.str();
}

// 辅助函数，检查文件是否存在
#include <sys/stat.h>
bool file_exists(const std::string& filename) {
    struct stat buffer;   
    return (stat(filename.c_str(), &buffer) == 0); 
}

// std::string generate_index_filename(
//     const std::string& method,
//     const std::string& index_dir,
//     size_t ref_len,
//     size_t anchor_len, // k
//     const MTree& config_tree // 传入一个配置好的 MTree
// ) {
//     std::stringstream ss;
//     ss << index_dir << "/" << method 
//        << "_ref" << ref_len 
//        << "_k" << anchor_len
//        << "_minN" << config_tree.minNodeCapacity
//        << "_maxN" << config_tree.maxNodeCapacity
//        << "_leafR" << config_tree.leafRadiusThreshold
//        << "_compMin" << config_tree.minCompactnessThreshold
//        << "_compFac" << config_tree.preSplitRadiusRatio;
        
//     ss << ".mtree_index";
//     return ss.str();
// }

// std::string get_fm_index_cache_path(size_t ref_len) {
//     // 直接使用传入的长度参数
//     return "fm_index_cache_" + std::to_string(ref_len) + ".bin";
// }

// // 写入缓存的函数
// bool write_fm_index_cache(const seqan::Index<seqan::Dna5String, seqan::FMIndex<>>& fm_index, size_t ref_len) {
//     const std::string cache_path = get_fm_index_cache_path(ref_len);

//     std::fstream file(cache_path, 
//                       std::ios::out | std::ios::binary);
//     if (!file.is_open()) {
//         std::cerr << "[ERROR] Failed to open cache file for writing: " << cache_path << "\n";
//         return false;
//     }
    
//     seqan::write(file, fm_index);
//     file.close();
//     std::cout << "[INFO] FM-Index cached successfully to " << cache_path << "\n";
//     return true;
// }

// // 读取缓存的函数
// bool load_fm_index_cache(seqan::Index<seqan::Dna5String, seqan::FMIndex<>>& fm_index, size_t ref_len) {
//     const std::string cache_path = get_fm_index_cache_path(ref_len);

//     std::fstream file(cache_path, 
//                       std::ios::in | std::ios::binary);
//     if (!file.is_open()) {
//         return false; // 文件不存在，表示没有缓存
//     }

//     try {
//         seqan::load(file, fm_index);
//         file.close();
//         std::cout << "[INFO] FM-Index loaded successfully from cache: " << cache_path << "\n";
//         return true;
//     } catch (const seqan::Exception& e) {
//         std::cerr << "[ERROR] Failed to load FM-Index from cache: " << e.what() << ". Rebuilding index.\n";
//         file.close();
//         // 如果加载失败，删除损坏的文件
//         std::remove(cache_path.c_str()); 
//         return false;
//     }
// }