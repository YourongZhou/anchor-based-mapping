#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_


#include <ext/algorithm>
#include <iterator>
#include <set>
#include <utility>
#include <vector>
#include <cmath>
#include <map>
#include <limits>
#include <unordered_map>
#include <algorithm>
#include <cassert>

#include "mtree_types.h"


namespace mt {
namespace functions {


/**
 * @brief A distance function object which calculates the <b>euclidean
 * distance</b> between two data objects representing coordinates.
 * @details Assumes that the data objects are same-sized sequences of numbers.
 * @see http://en.wikipedia.org/wiki/Euclidean_distance
 */
struct euclidean_distance {

	/**
	 * @brief  The operator that performs the calculation.
	 */
	template <typename Sequence>
	double operator()(const Sequence& data1, const Sequence& data2) const {
		double distance = 0;
		for(auto i1 = data1.begin(), i2 = data2.begin(); i1 != data1.end()  &&  i2 != data2.end(); ++i1, ++i2) {
			double diff = *i1 - *i2;
			distance += diff * diff;
		}
		distance = sqrt(distance);
		return distance;
	}
};



/**
 * @brief A promotion function object which randomly chooses two data objects
 * as promoted.
 */
struct random_promotion {
	/**
	 * @brief  The operator that performs the promotion.
	 * @tparam Data The type of the data objects.
	 * @tparam DistanceFunction The type of the function or function object used
	 *         to calculate the distance between two Data objects.
	 * @return A pair with the promoted data objects.
	 */
	template <typename Data, typename DistanceFunction>
	std::pair<Data, Data> operator()(const std::set<Data>& data_objects, DistanceFunction& distance_function) const {
		std::vector<Data> promoted;
		__gnu_cxx::random_sample_n(data_objects.begin(), data_objects.end(), inserter(promoted, promoted.begin()), 2);
		assert(promoted.size() == 2);
		return {promoted[0], promoted[1]};
	}
};



/**
 * @brief A partition function object which equally distributes the data objects
 *        according to their distances to the promoted data objects.
 * @details The algorithm is roughly equivalent to this:
 * @code
 *     data_objects := first_partition
 *     first_partition  := Empty
 *     second_partition := Empty
 *     Repeat until data_object is empty:
 *         X := The object in data_objects which is the nearest to promoted.first
 *         Remove X from data_object
 *         Add X to first_partition
 *
 *         Y := The object in data_objects which is the nearest to promoted.second
 *         Remove Y from data_object
 *         Add Y to second_partition
 * @endcode
 */
struct balanced_partition {
	/**
	 * @brief  The operator that performs the partition.
	 * @tparam Data The type of the data objects.
	 * @tparam DistanceFunction The type of the function or function object used
	 *                          to calculate the distance between two @c Data
	 *                          objects.
	 * @param [in]     promoted        The promoted data objects.
	 * @param [in,out] first_partition Initially, is the set containing all the
	 *                                 objects that must be partitioned. After
	 *                                 the partitioning, contains the objects
	 *                                 related to the first promoted data object,
	 *                                 which is @c promoted.first.
	 * @param [out]   second_partition Initially, is an empty set. After the
	 *                                 partitioning, contains the objects related
	 *                                 to the second promoted data object, which
	 *                                 is @c promoted.second.
	 * @param [in]   distance_function The distance function or function object.
	 */
	template <typename Data, typename DistanceFunction>
	void operator()(const std::pair<Data, Data>& promoted,
	                std::set<Data>& first_partition,
	                std::set<Data>& second_partition,
	                DistanceFunction& distance_function
	            ) const
	{
		std::vector<Data> queue1(first_partition.begin(), first_partition.end());
		// Sort by distance to the first promoted data
		std::sort(queue1.begin(), queue1.end(),
			[&](const Data& data1, const Data& data2) {
				double distance1 = distance_function(data1, promoted.first);
				double distance2 = distance_function(data2, promoted.first);
				return distance1 < distance2;
			}
		);

		std::vector<Data> queue2(first_partition.begin(), first_partition.end());
		// Sort by distance to the second promoted data
		std::sort(queue2.begin(), queue2.end(),
			[&](const Data& data1, const Data& data2) {
				double distance1 = distance_function(data1, promoted.second);
				double distance2 = distance_function(data2, promoted.second);
				return distance1 < distance2;
			}
		);

		first_partition.clear();

		typename std::vector<Data>::iterator i1 = queue1.begin();
		typename std::vector<Data>::iterator i2 = queue2.begin();

		while(i1 != queue1.end()  ||  i2 != queue2.end()) {
			while(i1 != queue1.end()) {
				Data& data = *i1;
				++i1;
				if(second_partition.find(data) == second_partition.end()) {
					first_partition.insert(data);
					break;
				}
			}

			while(i2 != queue2.end()) {
				Data& data = *i2;
				++i2;
				if(first_partition.find(data) == first_partition.end()) {
					second_partition.insert(data);
					break;
				}
			}
		}
	}
};

// --- 新的通用分裂结果类型 ---
/**
 * @brief 表示多路分裂的结果。
 * vector中的每个pair包含一个PromotedData（分区中心）和一个Partition（该分区的所有数据）。
 */
// template <typename Data>
// typedef std::set<Data> Partition;
// using SplitResult = std::vector<std::pair<Data, Partition>>;

/**
 * @brief A function object that defines a split function by composing a
 *        promotion function and a partition function.
 * @tparam PromotionFunction The type of the function or function object which
 *                           implements a promotion function.
 * @tparam PartitionFunction The type of the function or function object which
 *                           implements a partition function.
 */
// template <typename PromotionFunction, typename PartitionFunction>
// struct split_function {
// 	/** */
// 	typedef PromotionFunction promotion_function_type;

// 	/** */
// 	typedef PartitionFunction partition_function_type;

// 	PromotionFunction promotion_function;
// 	PartitionFunction partition_function;

// 	/** */
// 	explicit split_function(
// 			PromotionFunction promotion_function = PromotionFunction(),
// 			PartitionFunction partition_function = PartitionFunction()
// 		)
// 	: promotion_function(promotion_function),
// 	  partition_function(partition_function)
// 	{}


// 	/**
// 	 * @brief The operator that performs the split.
// 	 * @tparam Data The type of the data objects.
// 	 * @tparam DistanceFunction The type of the function or function object used
// 	 *                          to calculate the distance between two @c Data
// 	 *                          objects.
// 	 * @param [in,out] first_partition Initially, is the set containing all the
// 	 *                                 objects that must be partitioned. After
// 	 *                                 the partitioning, contains the objects
// 	 *                                 related to the first promoted data object.
// 	 * @param [out]   second_partition Initially, is an empty set. After the
// 	 *                                 partitioning, contains the objects related
// 	 *                                 to the second promoted data object.
// 	 * @param [in]   distance_function The distance function or function object.
// 	 * @return A pair with the promoted data objects.
// 	 */
// 	template <typename Data, typename DistanceFunction>
// 	std::pair<Data, Data> operator()(
// 				std::set<Data>& first_partition,
// 				std::set<Data>& second_partition,
// 				DistanceFunction& distance_function
// 			) const
// 	{
// 		// --- promotion ---
// 		std::pair<Data, Data> promoted = promotion_function(first_partition, distance_function);

// 		// --- 记录原始的全部数据点，用于统计 ---
// 		std::vector<Data> all_data(first_partition.begin(), first_partition.end());

// 		// --- partition ---
// 		partition_function(promoted, first_partition, second_partition, distance_function);

// 		// --- 统计部分 ---
// 		size_t total_points = all_data.size();
// 		size_t first_count = first_partition.size();
// 		size_t second_count = second_partition.size();

// 		// 计算 range
// 		double range_first = 0.0;
// 		for (const auto& d : first_partition)
// 			range_first = std::max(range_first, distance_function(d, promoted.first));

// 		double range_second = 0.0;
// 		for (const auto& d : second_partition)
// 			range_second = std::max(range_second, distance_function(d, promoted.second));

// 		// 统计“落在另一边 range 内”的点
// 		size_t cross_to_second = 0;
// 		for (const auto& d : second_partition) {
// 			double dist1 = distance_function(d, promoted.first);
// 			if (dist1 <= range_first)
// 				++cross_to_second;
// 		}

// 		size_t cross_to_first = 0;
// 		for (const auto& d : first_partition) {
// 			double dist2 = distance_function(d, promoted.second);
// 			if (dist2 <= range_second)
// 				++cross_to_first;
// 		}

// 		std::cout << "---- Split stats ----\n";
// 		std::cout << "Total points: " << total_points << "\n";
// 		std::cout << "First partition: " << first_count << "\n";
// 		std::cout << "Second partition: " << second_count << "\n";
// 		std::cout << "Cross to second: " << cross_to_second << "\n";
// 		std::cout << "Cross to first: " << cross_to_first << "\n";
// 		std::cout << "---------------------\n";

// 		// --- 返回 ---
// 		return promoted;
// 	}
// };

// 定义一个内部 Cluster 结构体 (用于辅助聚类)
const double MAX_CAP_RATIO = 0.60;
const int MIN_CLUSTERS = 2; // 至少要分裂成两个
const double MAX_RADIUS_RATIO = 1.5;

// 定义距离矩阵类型
using DistanceMatrix = std::vector<std::vector<double>>;
template <typename DataT>
struct InternalCluster {
    int center_index;              // 推广数据的索引
    std::set<int> entry_indices;   // Cluster 中包含的所有数据点的索引
    double radius;                 // 覆盖半径
};

template <typename DataT>
void recalculate_center_and_radius(
    InternalCluster<DataT>& cluster,
    const DistanceMatrix& dist_matrix
) {
    if (cluster.entry_indices.empty()) {
        cluster.radius = 0.0;
        return;
    }

    double best_max_dist = std::numeric_limits<double>::max();
    int new_center_index = -1;

    // 遍历所有可能的中心点 (即 Cluster 中的每个条目)
    for (const auto& candidate_center_idx : cluster.entry_indices) {
        double current_max_dist = 0.0;

        // 找到该中心到 Cluster 中最远点的距离 (即覆盖半径)
        for (const auto& entry_idx : cluster.entry_indices) {
            // !!! 核心优化：查表 !!!
            double dist = dist_matrix[candidate_center_idx][entry_idx];
            
            if (dist > current_max_dist) {
                current_max_dist = dist;
            }
        }

        // M-Tree 启发式：选择使覆盖半径最小的点作为中心
        if (current_max_dist < best_max_dist) {
            best_max_dist = current_max_dist;
            new_center_index = candidate_center_idx;
        }
    }

    cluster.radius = best_max_dist;
    if (new_center_index != -1) {
        cluster.center_index = new_center_index;
    }
}

// optimized_split 函数
template <typename Data, typename DistanceFunction>
SplitResult<Data> optimized_split(
    std::set<Data>& all_entries,
    DistanceFunction& distance_function,
    double current_radius // 分裂前节点的原始半径 (R_orig)
) {
    const size_t TOTAL_SIZE = all_entries.size();
    if (TOTAL_SIZE == 0) return {};

    // 1. 距离矩阵预计算 (空间换时间)

    // A. 建立 Data <-> Index 的映射
    std::unordered_map<Data, int> data_to_index;
    std::vector<Data> index_to_data;
    int index = 0;
    for (const auto& entry : all_entries) {
        data_to_index[entry] = index;
        index_to_data.push_back(entry);
        index++;
    }
    const size_t N = TOTAL_SIZE;

    // B. 构建距离矩阵
    DistanceMatrix dist_matrix(N, std::vector<double>(N));
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i; j < N; ++j) {
            double dist = distance_function(index_to_data[i], index_to_data[j]);
            dist_matrix[i][j] = dist;
            dist_matrix[j][i] = dist; // 矩阵对称
        }
    }

    // 2. 初始化：创建原子 Cluster 集合
    std::vector<InternalCluster<Data>> clusters;
    clusters.reserve(N);

    for (int i = 0; i < N; ++i) {
        InternalCluster<Data> c;
        c.center_index = i; // 初始中心就是自己
        c.entry_indices.insert(i);
        c.radius = 0.0;
        clusters.push_back(std::move(c));
    }

    // 3. 贪婪合并循环
    const size_t MAX_CLUSTER_SIZE = static_cast<size_t>(TOTAL_SIZE * MAX_CAP_RATIO);
    const double MAX_ALLOWED_RADIUS = current_radius * MAX_RADIUS_RATIO;
    
    // 循环条件：Cluster 数量大于最小要求 (MIN_CLUSTERS)
    while (clusters.size() > MIN_CLUSTERS) {
        double min_dist = std::numeric_limits<double>::max();
        int best_idx1 = -1, best_idx2 = -1;

        // 3a. 寻找最佳合并对 (最近的 Cluster 中心)
        for (int i = 0; i < clusters.size(); ++i) {
            for (int j = i + 1; j < clusters.size(); ++j) {
                // !!! 核心优化：查表 !!!
                double dist = dist_matrix[clusters[i].center_index][clusters[j].center_index];
                if (dist < min_dist) {
                    min_dist = dist;
                    best_idx1 = i;
                    best_idx2 = j;
                }
            }
        }
        
        if (best_idx1 == -1) break; // 找不到任何一对

        // 3b. 计算试探性 Cluster C_new
        InternalCluster<Data> c_new;
        
        // 合并条目索引
        c_new.entry_indices.insert(clusters[best_idx1].entry_indices.begin(), clusters[best_idx1].entry_indices.end());
        c_new.entry_indices.insert(clusters[best_idx2].entry_indices.begin(), clusters[best_idx2].entry_indices.end());
        
        // 重新计算中心和半径 (使用查表)
        recalculate_center_and_radius(c_new, dist_matrix);

        // 3c. 约束检查
        bool capacity_ok = c_new.entry_indices.size() <= MAX_CLUSTER_SIZE;
        bool radius_ok = c_new.radius <= MAX_ALLOWED_RADIUS;

        // 3d. 执行合并或终止
        if (capacity_ok && radius_ok) {
            // 合并：移除旧的，添加新的
            
            // 确保移除索引大的 (best_idx2) 优先，防止索引失效
            clusters.erase(clusters.begin() + best_idx2); 
            clusters.erase(clusters.begin() + best_idx1);
            
            clusters.push_back(std::move(c_new));
            
        } else {
            // 终止合并循环
            break; 
        }
    }
    
    // 4. 转换为 SplitResult 并清空输入

    SplitResult<Data> result;
    result.reserve(clusters.size());
    
    // 遍历最终的 Cluster 结果
    for (const auto& cluster : clusters) {
        Partition final_partition;
        
        // 将索引转换回 Data 对象
        for (const auto& idx : cluster.entry_indices) {
            final_partition.insert(index_to_data[idx]);
        }
        
        // 确保中心点的数据对象是正确的
        const Data& promoted_data = index_to_data[cluster.center_index];
        
        result.push_back({promoted_data, final_partition});
    }

    all_entries.clear(); // 清空原始输入集合
    
    return result;
}

/**
 * @brief 策略实现：使用外部优化的多路分裂函数。
 */
struct OptimizedKSplitStrategy {

    template <typename Data, typename DistanceFunction>
    SplitResult<Data> operator()(
        std::set<Data>& all_entries,
        DistanceFunction& distance_function,
		double current_radius
    ) const {
        return optimized_split(all_entries, distance_function, current_radius);
    }
};

/**
 * @brief 兼容策略实现：将现有的二路分裂包装成 Multi-Split 签名。
 */
template <typename PromotionFunction, typename PartitionFunction>
struct TwoWaySplitStrategy {
    
    PromotionFunction promotion_function;
    PartitionFunction partition_function;

    explicit TwoWaySplitStrategy(
            PromotionFunction pf = PromotionFunction(),
            PartitionFunction paf = PartitionFunction()
        )
    : promotion_function(pf), partition_function(paf)
    {}

    template <typename DataT, typename DistanceFunction>
    SplitResult<DataT> operator()(
        std::set<DataT>& all_entries,
        DistanceFunction& distance_function,
		double current_radius
    ) const {
        // 1. 复制一份 all_entries 作为 first_partition 的初始值
        //    (因为原来的签名是 in/out 参数)
        Partition first_partition = all_entries;
        Partition second_partition;

        // 2. 调用现有的 promotion 和 partition 逻辑
        std::pair<Data, Data> promoted = promotion_function(first_partition, distance_function);

		// --- 记录原始的全部数据点，用于统计 ---
		std::vector<Data> all_data(first_partition.begin(), first_partition.end());

        partition_function(promoted, first_partition, second_partition, distance_function);

        // 3. 将结果转换成新的多路格式
        SplitResult<Data> result;
        result.push_back({promoted.first, first_partition});
        result.push_back({promoted.second, second_partition});

        // 4. 清空原始的 all_entries，因为数据已经被分配到 result 中
        all_entries.clear(); 

		// --- 统计部分 ---
		size_t total_points = all_data.size();
		size_t first_count = first_partition.size();
		size_t second_count = second_partition.size();

		// 计算 range
		double range_first = 0.0;
		for (const auto& d : first_partition)
			range_first = std::max(range_first, distance_function(d, promoted.first));

		double range_second = 0.0;
		for (const auto& d : second_partition)
			range_second = std::max(range_second, distance_function(d, promoted.second));

		// 统计“落在另一边 range 内”的点
		size_t cross_to_second = 0;
		for (const auto& d : second_partition) {
			double dist1 = distance_function(d, promoted.first);
			if (dist1 <= range_first)
				++cross_to_second;
		}

		size_t cross_to_first = 0;
		for (const auto& d : first_partition) {
			double dist2 = distance_function(d, promoted.second);
			if (dist2 <= range_second)
				++cross_to_first;
		}

		std::cout << "---- Split stats ----\n";
		std::cout << "Total points: " << total_points << "\n";
		std::cout << "First partition: " << first_count << "\n";
		std::cout << "Second partition: " << second_count << "\n";
		std::cout << "Cross to second: " << cross_to_second << "\n";
		std::cout << "Cross to first: " << cross_to_first << "\n";
		std::cout << "---------------------\n";
        
        return result;
    }
};

template <typename Data, typename DistanceFunction>
class cached_distance_function {
public:
	explicit cached_distance_function(const DistanceFunction& distance_function)
		: distance_function(distance_function)
		{}

	double operator()(const Data& data1, const Data& data2) {
		typename CacheType::iterator i = cache.find(std::make_pair(data1, data2));
		if(i != cache.end()) {
			return i->second;
		}

		i = cache.find(std::make_pair(data2, data1));
		if(i != cache.end()) {
			return i->second;
		}

		// Not found in cache
		double distance = distance_function(data1, data2);

		// Store in cache
		cache.insert(std::make_pair(std::make_pair(data1, data2), distance));
		cache.insert(std::make_pair(std::make_pair(data2, data1), distance));

		return distance;
	}

private:
	typedef std::map<std::pair<Data, Data>, double> CacheType;

	const DistanceFunction& distance_function;
	CacheType cache;
};

template <typename Data>
struct DataPtrLess {
    // This assumes your Data type has a defined operator<
    bool operator()(const Data* lhs, const Data* rhs) const {
        // Check for nullptr safety if needed, but for M-Tree pivots, they 
        // shouldn't be null here.
        return *lhs < *rhs;
    }
};



} /* namespace functions */
} /* namespace mtree */


#endif /* FUNCTIONS_H_ */
