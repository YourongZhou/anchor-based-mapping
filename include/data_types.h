// ----------------- include/data_types.h -----------------
#ifndef DATA_TYPES_H
#define DATA_TYPES_H

using namespace std;
using namespace mt;

// Levenshtein 距离封装，只比较 seq
struct SubstringLevDist {
    double operator()(const Substring &a, const Substring &b) const {
        return static_cast<double>(levenshtein(a.seq, b.seq));
    }
};

using Data = Substring;
using Distance = SubstringLevDist;
using CachedDistance = functions::cached_distance_function<Data, Distance>;

// split function 类型定义
// using SplitStrategyType = functions::TwoWaySplitStrategy<
//     functions::random_promotion,
//     functions::balanced_partition
// >;
using SplitStrategyType = functions::OptimizedKSplitStrategy;
using MTree = mtree<Data, Distance, SplitStrategyType>;

struct CandidateResults {
    std::vector<int> positions; // 原函数返回的起始位置
    double average_distance;    // 平均编辑距离
};

#endif // DATA_TYPES_H