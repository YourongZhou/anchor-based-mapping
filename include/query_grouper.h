#ifndef QUERY_GROUPER_H
#define QUERY_GROUPER_H

#include <vector>
#include <string>
#include <limits>
#include <omp.h>
#include "data_types.h"

/**
 * @brief QueryGrouper: 查询分桶器，用于 Cache-Aware WORM 策略
 * @details 将查询按几何距离分组到不同的锚点桶中，
 *          同一桶内的查询会访问 M-Tree 的同一子树，提高缓存复用率
 */
class QueryGrouper {
public:
    /**
     * @brief 分组批次结构
     */
    struct GroupedBatch {
        int anchor_index;                    // 锚点索引
        std::vector<int> query_indices;      // 属于该桶的查询索引（存索引而非字符串，节省内存）
    };

    /**
     * @brief 核心分桶逻辑：将查询分配到最近的锚点桶
     * @param queries 查询字符串向量
     * @param anchors 路由锚点向量（来自 M-Tree 顶层节点）
     * @param distFunc 距离函数对象
     * @return 分组后的批次向量，每个批次包含一个锚点索引和属于该桶的查询索引列表
     */
    static std::vector<GroupedBatch> groupQueries(
        const std::vector<std::string>& queries, 
        const std::vector<Data>& anchors,
        Distance& distFunc
    ) {
        if (anchors.empty() || queries.empty()) {
            return {};
        }

        int num_anchors = static_cast<int>(anchors.size());
        std::vector<GroupedBatch> batches(num_anchors);
        for(int i = 0; i < num_anchors; ++i) {
            batches[i].anchor_index = i;
        }

        // 使用 OpenMP 并行预处理，因为是纯计算，不涉及 IO
        // 每个线程维护私有的批次向量，最后合并
        #pragma omp parallel
        {
            // 每个线程的私有批次向量
            std::vector<std::vector<int>> thread_batches(num_anchors);
            for(int i = 0; i < num_anchors; ++i) {
                thread_batches[i].reserve(queries.size() / (num_anchors * omp_get_num_threads()) + 1);
            }

            #pragma omp for nowait
            for (size_t i = 0; i < queries.size(); ++i) {
                double min_dist = std::numeric_limits<double>::max();
                int best_anchor = 0;

                // 将查询字符串转换为 Data 类型
                Data query_data{queries[i], -1}; // pos = -1 表示查询
                
                // 寻找最近的 Top-Level Anchor
                for (int j = 0; j < num_anchors; ++j) {
                    double d = distFunc(query_data, anchors[j]);
                    if (d < min_dist) {
                        min_dist = d;
                        best_anchor = j;
                    }
                }
                
                // 将查询索引添加到对应线程的批次中
                thread_batches[best_anchor].push_back(static_cast<int>(i));
            }

            // 合并线程私有批次到全局批次（临界区）
            #pragma omp critical
            {
                for (int i = 0; i < num_anchors; ++i) {
                    batches[i].query_indices.insert(
                        batches[i].query_indices.end(),
                        thread_batches[i].begin(),
                        thread_batches[i].end()
                    );
                }
            }
        }

        // 移除空的批次
        std::vector<GroupedBatch> non_empty_batches;
        non_empty_batches.reserve(batches.size());
        for (auto& batch : batches) {
            if (!batch.query_indices.empty()) {
                non_empty_batches.push_back(std::move(batch));
            }
        }

        return non_empty_batches;
    }
};

#endif // QUERY_GROUPER_H

