#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include <vector>
#include <string>
#include "mtree_types.h"
#include "functions.h"
#include "levenshtein.hpp"
#include "mtree.h"  // 必须包含 mtree.h 才能定义 MTree 别名

// Substring 已经在 mtree_types.h 中定义

using namespace mt;

// Levenshtein 距离封装，只比较 seq
struct SubstringLevDist {
    double operator()(const Substring &a, const Substring &b) const {
        return static_cast<double>(levenshtein(a.seq, b.seq));
    }
};

using Data = Substring;
using Distance = SubstringLevDist;
using CachedDistance = mt::functions::cached_distance_function<Data, Distance>;

// split function 类型定义
using SplitStrategyType = mt::functions::OptimizedKSplitStrategy;
using MTree = mt::mtree<Data, Distance, SplitStrategyType>;

struct CandidateResults {
    std::vector<int> positions; // 原函数返回的起始位置
    double average_distance;    // 平均编辑距离
    // 增加统计信息字段
    size_t node_access = 0;
    std::vector<double> leaf_node_radii;
};

#endif // DATA_TYPES_H
