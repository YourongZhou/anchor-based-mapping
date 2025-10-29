// ----------------- include/mtree_types.h -----------------
#ifndef MTREE_TYPES_H
#define MTREE_TYPES_H

#include <set>
#include <vector>
#include <utility>
#include <functional>


// 假设 Substring 是一个已定义的类
// ======================== 类型定义 ========================
struct Substring {
    std::string seq;
    int pos;

    bool operator==(const Substring &other) const {
        return seq == other.seq && pos == other.pos;
    }

    bool operator<(const Substring& other) const {
        return seq < other.seq;
    }
};

namespace std {
    // 特化 std::hash 模板，使其适用于 Substring 类型
    template <>
    struct hash<Substring> {
        size_t operator()(const Substring& s) const noexcept {
            // 1. 获取 std::string 的哈希对象
            std::hash<std::string> string_hasher;
            
            // 2. 获取 int 的哈希对象
            std::hash<int> int_hasher;
            
            // 3. 计算各个部分的哈希值
            size_t h1 = string_hasher(s.seq);
            size_t h2 = int_hasher(s.pos);

            // 4. 组合哈希值 (使用 boost::hash_combine 的思想)
            // 这是一个标准的哈希组合公式：h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2))
            // C++20 提供了 std::hash_combine，但使用一个简化的版本更兼容：
            return h1 ^ (h2 << 1); // 简单且兼容的组合方式

            /* // 推荐的更健壮的组合方式 (类似 boost::hash_combine):
            // 更好的组合哈希值的函数，可以放在一个工具函数中
            size_t seed = h1;
            seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
            */
        }
    };
}

// --- 基础类型定义 ---
typedef Substring Data; 
typedef std::set<Data> Partition; // <--- 从 mtree.h 移到这里

// --- 通用分裂结果类型定义（依赖于 Data 和 Partition）---
template <typename T>
using SplitResult = std::vector<std::pair<T, Partition>>;

#endif // MTREE_TYPES_H