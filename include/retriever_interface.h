#ifndef RETRIEVER_INTERFACE_H
#define RETRIEVER_INTERFACE_H

#include <string>
#include <vector>
#include <memory>
#include "data_types.h"

// 抽象检索器接口
class ICandidateRetriever {
public:
    virtual ~ICandidateRetriever() = default;

    // 构建索引
    virtual void buildIndex(const std::string& reference) = 0;

    // 检索候选位置
    // query: 查询序列
    // maxDist: 最大允许编辑距离
    // requireDistance: 是否需要计算平均编辑距离
    virtual CandidateResults retrieve(const std::string& query, int maxDist, bool requireDistance) = 0;

    // 报告算法特有的统计信息
    virtual void reportStats() const {}
};

#endif // RETRIEVER_INTERFACE_H

