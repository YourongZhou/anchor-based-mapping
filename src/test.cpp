#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include "mtree.h"
#include "functions.h"

// 定义别名以简化代码
using namespace mt;
using std::vector;
using std::cout;
using std::endl;

// 数据类型
using Data = vector<double>;

// 定义一个 MTree 类型
using Distance = functions::euclidean_distance;
using SplitFunc = functions::split_function<
    std::pair<Data, Data>(*)(const std::set<Data>&, functions::cached_distance_function<Data, Distance>&),
    functions::balanced_partition>;
using MTree = mtree<Data, Distance, SplitFunc>;

// 简单的非随机 promotion 函数
std::pair<Data, Data> nonRandomPromotion(
    const std::set<Data>& dataSet,
    functions::cached_distance_function<Data, Distance>&)
{
    vector<Data> objs(dataSet.begin(), dataSet.end());
    std::sort(objs.begin(), objs.end());
    return {objs.front(), objs.back()};
}

int main() {
    cout << "===== MTree Demo =====" << endl;

    // 创建 MTree
    MTree tree(
        2,                // min node capacity
        -1,               // max node capacity (-1 表示自动推导)
        Distance(),       // 距离函数
        SplitFunc(nonRandomPromotion) // split 函数
    );

    // 添加数据点
    tree.add({0.0, 0.0});
    tree.add({1.0, 1.0});
    tree.add({2.0, 2.0});
    tree.add({5.0, 5.0});
    tree.add({10.0, 10.0});

    cout << "Inserted 5 points into the tree." << endl;

    // 查询最近邻（半径内搜索）
    Data query = {1.5, 1.5};
    double radius = 3.0;
    auto results = tree.get_nearest_by_range(query, radius);

    cout << "Query point: (" << query[0] << ", " << query[1] << ")\n";
    cout << "Searching within radius " << radius << "...\n";

    for (auto it = results.begin(); it != results.end(); ++it) {
        cout << "  Found: (" << it->data[0] << ", " << it->data[1]
             << ")  dist=" << it->distance << endl;
    }

    // 测试删除
    tree.remove({2.0, 2.0});
    cout << "Removed (2, 2) from the tree.\n";

    // 再次搜索
    results = tree.get_nearest_by_limit(query, 3);
    cout << "Top 3 nearest points now:\n";
    for (auto it = results.begin(); it != results.end(); ++it) {
        cout << "  (" << it->data[0] << ", " << it->data[1]
             << ")  dist=" << it->distance << endl;
    }

    cout << "===== DONE =====" << endl;
    return 0;
}
