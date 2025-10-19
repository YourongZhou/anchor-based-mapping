#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_


#include <ext/algorithm>
#include <iterator>
#include <set>
#include <utility>
#include <vector>
#include <cmath>
#include <map>


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



/**
 * @brief A function object that defines a split function by composing a
 *        promotion function and a partition function.
 * @tparam PromotionFunction The type of the function or function object which
 *                           implements a promotion function.
 * @tparam PartitionFunction The type of the function or function object which
 *                           implements a partition function.
 */
template <typename PromotionFunction, typename PartitionFunction>
struct split_function {
	/** */
	typedef PromotionFunction promotion_function_type;

	/** */
	typedef PartitionFunction partition_function_type;

	PromotionFunction promotion_function;
	PartitionFunction partition_function;

	/** */
	explicit split_function(
			PromotionFunction promotion_function = PromotionFunction(),
			PartitionFunction partition_function = PartitionFunction()
		)
	: promotion_function(promotion_function),
	  partition_function(partition_function)
	{}


	/**
	 * @brief The operator that performs the split.
	 * @tparam Data The type of the data objects.
	 * @tparam DistanceFunction The type of the function or function object used
	 *                          to calculate the distance between two @c Data
	 *                          objects.
	 * @param [in,out] first_partition Initially, is the set containing all the
	 *                                 objects that must be partitioned. After
	 *                                 the partitioning, contains the objects
	 *                                 related to the first promoted data object.
	 * @param [out]   second_partition Initially, is an empty set. After the
	 *                                 partitioning, contains the objects related
	 *                                 to the second promoted data object.
	 * @param [in]   distance_function The distance function or function object.
	 * @return A pair with the promoted data objects.
	 */
	template <typename Data, typename DistanceFunction>
	std::pair<Data, Data> operator()(
				std::set<Data>& first_partition,
				std::set<Data>& second_partition,
				DistanceFunction& distance_function
			) const
	{
		// --- promotion ---
		std::pair<Data, Data> promoted = promotion_function(first_partition, distance_function);

		// --- 记录原始的全部数据点，用于统计 ---
		std::vector<Data> all_data(first_partition.begin(), first_partition.end());

		// --- partition ---
		partition_function(promoted, first_partition, second_partition, distance_function);

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

		// --- 返回 ---
		return promoted;
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



} /* namespace functions */
} /* namespace mtree */


#endif /* FUNCTIONS_H_ */
