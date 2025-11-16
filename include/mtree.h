#ifndef MTREE_H_
#define MTREE_H_

#include <iterator>
#include <limits>
#include <map>
#include <queue>
#include <utility>
#include <cassert>

#include "mtree_types.h"
#include "functions.h"
#include "access_log.hpp"



namespace mt {



/**
 * @brief The main class that implements the M-Tree.
 *
 * @tparam Data The type of data that will be indexed by the M-Tree. This type
 *         must be an assignable type and a strict weak ordering must be defined
 *         by @c std::less<Data>.
 * @tparam DistanceFunction The type of the function that will be used to
 *         calculate the distance between two @c Data objects. By default, it is
 *         ::mt::functions::euclidean_distance.
 * @tparam SplitFunction The type of the function that will be used to split a
 *         node when it is at its maximum capacity and a new child must be
 *         added. By default, it is a composition of
 *         ::mt::functions::random_promotion and
 *         ::mt::functions::balanced_partition.
 *
 *
 * @todo Include a @c Compare template and constructor parameters instead of
 *       implicitly using @c std::less<Data> on @c stc::set and @c std::map.
 *
 * @todo Implement an @c unordered_mtree class which uses @c std::unordered_set
 *      and @c std::unordered_map instead of @c std::set and @c std::map
 *      respectively.
 */
template <
	typename Data,
	typename DistanceFunction = ::mt::functions::euclidean_distance,
	typename SplitStrategy = ::mt::functions::TwoWaySplitStrategy<
	        ::mt::functions::random_promotion,
	        ::mt::functions::balanced_partition
		>
	// typename SplitFunction = ::mt::functions::split_function<
	//         ::mt::functions::random_promotion,
	//         ::mt::functions::balanced_partition
	// 	>
>
class mtree {
public:
	typedef DistanceFunction distance_function_type;
	typedef SplitStrategy split_strategy_type;
	// typedef SplitFunction    split_function_type;
	typedef functions::cached_distance_function<Data, DistanceFunction> cached_distance_function_type;
	// --- 紧凑性分裂的配置参数 ---
    size_t minCompactnessThreshold; // 触发紧凑性分裂的最小节点容量
    double preSplitRadiusRatio;     // 触发紧凑性分裂的半径膨胀比例

private:
	class Node;
	class Entry;

	// Exception classes
	class SplitNodeReplacement {
	public:
		// 使用 std::vector 存储新节点，支持任意数量
		std::vector<Node*> newNodes;

		/**
		 * @brief 构造函数：接受一个包含所有新节点的 vector。
		 * @param newNodes 包含所有 k 个新节点的 vector。
		 */
		SplitNodeReplacement(std::vector<Node*> newNodes) 
			// 使用 std::move 提高效率，避免不必要的深拷贝
			: newNodes(std::move(newNodes)) 
		{
			// 确保至少有两个新节点（因为是分裂操作）
			if (this->newNodes.size() < 2) {
				std::cerr << "Error: SplitNodeReplacement created with less than 2 nodes." << std::endl;
			}
		}
		
		// 获取新节点数量
		size_t numNewNodes() const {
			return newNodes.size();
		}
	};


	// class SplitNodeReplacement {
	// public:
	// 	enum { NUM_NODES = 2 };
	// 	Node* newNodes[NUM_NODES];
		// SplitNodeReplacement(Node* newNodes[NUM_NODES]) {
		// 	for(int i = 0; i < NUM_NODES; ++i) {
		// 		this->newNodes[i] = newNodes[i];
		// 	}

		// }
	// };

	class CompactnessPromotion {
	public:
		Node* promotedNode; // 携带被提升的新节点
		explicit CompactnessPromotion(Node* new_node) : promotedNode(new_node) {}
	};

	class RootNodeReplacement {
	public:
		Node* newRoot;
	};

	class NodeUnderCapacity { };

	class DataNotFound {
	public:
		Data data;
	};


public:

	/**
	 * @brief A container-like class which can be iterated to fetch the results
	 *        of a nearest-neighbors query.
	 * @details The neighbors are presented in non-decreasing order from the
	 *          @c query_data argument to the mtree::get_nearest() call.
	 *
	 *          The query on the M-Tree is executed during the iteration, as the
	 *          results are fetched. It means that, by the time when the @a n-th
	 *          result is fetched, the next result may still not be known, and
	 *          the resources allocated were only the necessary to identify the
	 *          @a n first results.
	 *
	 *          The objects in the container are mtree::query_result instances,
	 *          which contain a data object and the distance from the query
	 *          data object.
	 * @see mtree::get_nearest()
	 */
	class query {
	public:
		mutable size_t nodeAccess;
		mutable std::vector<double> accessedLeafNodeRadii;

		/**
		 * @brief The type of the results for nearest-neighbor queries.
		 */
		class result_item {
		public:
			/** @brief A nearest-neighbor */
			Data data;

			/** @brief The distance from the nearest-neighbor to the query data
			 *         object parameter.
			 */
			double distance;

			/** @brief Default constructor */
			result_item() = default;

			/** @brief Copy constructor */
			result_item(const result_item&) = default;

			/** @brief Move constructor */
			result_item(result_item&&) = default;

			~result_item() = default;

			/** @brief Copy assignment */
			result_item& operator=(const result_item&) = default;

			/** @brief Move assignment */
			result_item& operator=(result_item&& ri) {
				if(this != &ri) {
					data = std::move(ri.data);
					distance = ri.distance;
				}
				return *this;
			}
		};


		typedef result_item value_type;


		query() = delete;

		/**
		 * @brief Copy constructor.
		 */
		query(const query&) = default;

		/**
		 * @brief Move constructor.
		 */
		query(query&&) = default;

		query(const mtree* _mtree, const Data& data, double range, size_t limit)
			: _mtree(_mtree), data(data), range(range), limit(limit), nodeAccess(0)
			{}


		/** @brief Copy assignment. */
		query& operator=(const query&) = default;

		/** @brief Move assignment. */
		query& operator=(query&& q) {
			if(this != &q) {
				this->_mtree = q._mtree;
				this->range = q.range;
				this->limit = q.limit;
				this->data = std::move(q.data);
			}
			return *this;
		}



		/**
		 * @brief The iterator for accessing the results of nearest-neighbor
		 *        queries.
		 */
		class iterator {
		public:
			typedef std::input_iterator_tag iterator_category;
			typedef result_item             value_type;
			typedef signed long int         difference_type;
			typedef result_item*            pointer;
			typedef result_item&            reference;


			iterator() : isEnd(true) {}


			explicit iterator(const query* _query)
				: _query(_query),
				  isEnd(false),
				  yieldedCount(0)
			{
				_query->nodeAccess = 0; //每次开始迭代时重复计数
				_query->accessedLeafNodeRadii.clear();
				if(_query->_mtree->root == NULL) {
					isEnd = true;
					return;
				}

				double distance = _query->_mtree->distance_function(_query->data, _query->_mtree->root->data);
				double minDistance = std::max(distance - _query->_mtree->root->radius, 0.0);

				pendingQueue.push({_query->_mtree->root, distance, minDistance});
				nextPendingMinDistance = minDistance;

				fetchNext();
			}


			/** @brief Copy constructor. */
			iterator(const iterator&) = default;

			/** @brief Move constructor. */
			iterator(iterator&&) = default;

			~iterator() = default;

			/** @brief Copy assignment. */
			iterator& operator=(const iterator&) = default;

			/** @brief Move assignment. */
			iterator& operator=(iterator&& i) {
				if(this != &i) {
					this->_query = i._query;
					this->currentResultItem = std::move(i.currentResultItem);
					this->isEnd = i.isEnd;
					this->pendingQueue = std::move(i.pendingQueue);
					this->nextPendingMinDistance = i.nextPendingMinDistance;
					this->nearestQueue = std::move(i.nearestQueue);
					this->yieldedCount = i.yieldedCount;
				}
				return *this;
			}


			bool operator==(const iterator& ri) const {
				if(this->isEnd  &&  ri.isEnd) {
					return true;
				}

				if(this->isEnd  ||  ri.isEnd) {
					return false;
				}

				return  this->_query == ri._query
				    &&  this->yieldedCount == ri.yieldedCount;
			}

			bool operator!=(const iterator& ri) const {
				return ! this->operator==(ri);
			}


			/**
			 * @brief Advance the iterator to the next result.
			 */
			//@{
			// prefix
			iterator& operator++() {
				fetchNext();
				return *this;
			}

			// postfix
			iterator operator++(int) {
				iterator aCopy = *this;
				operator++();
				return aCopy;
			}
			//@}


			/**
			 * @brief Gives access to the current query result.
			 * @details An iterator instance always refers to the same
			 *          result_item instance; only the fields of the result_item
			 *          are updated when the iterator changes.
			 */
			//@{
			const result_item& operator*() const {
				return currentResultItem;
			}

			const result_item* operator->() const {
				return &currentResultItem;
			}
			//@}

		private:
			template <typename U>
			struct ItemWithDistances {
				const U* item;
				double distance;
				double minDistance;

				ItemWithDistances(const U* item, double distance, double minDistance)
					: item(item), distance(distance), minDistance(minDistance)
					{ }

				bool operator<(const ItemWithDistances& that) const {
					return (this->minDistance > that.minDistance);
				}

			};

			void fetchNext() {
				assert(! isEnd);

				if(isEnd  ||  yieldedCount >= _query->limit) {
					isEnd = true;
					return;
				}

				 while(!pendingQueue.empty()  ||  !nearestQueue.empty()) {
					if(prepareNextNearest()) {
						return;
					}

					assert(!pendingQueue.empty());

					ItemWithDistances<Node> pending = pendingQueue.top();
					pendingQueue.pop();
					// _query->nodeAccess++; // 计数 nodeAccess

					const Node* node = pending.item;

					// 遍历父节点的子节点（可能是 LeafNode 或 InternalNode）
                    for(auto i = node->children.begin(); i != node->children.end(); ++i) {
                        
                        // 1. 获取子节点指针
                        // 假设 ChildrenMap 是 <Data, std::shared_ptr<IndexItem>>
                        IndexItem* child = i->second; 

                        // 2. 第一次剪枝检查 (基于父节点)
                        if(std::abs(pending.distance - child->distanceToParent) - child->radius <= _query->range) {
                            
                            // 3. 计算昂贵的距离 (只对 Node/LeafNode 计算)
                            double childDistance = _query->_mtree->distance_function(_query->data, child->data);
                            double childMinDistance = std::max(childDistance - child->radius, 0.0);
                            
                            // 4. 第二次剪枝检查 (基于子节点)
                            if(childMinDistance <= _query->range) {
                                
                                LeafNode* leafnode = dynamic_cast<LeafNode*>(child);
                                
                                if(leafnode != NULL) {
                                    // 找到了一个相关的 LeafNode (L)
                                    _query->nodeAccess++; 
                                    _query->accessedLeafNodeRadii.push_back(leafnode->radius);

                                    // d(Anchor, Query) (即 d(P_L, Q))
                                    double d_anchor_query = childDistance;

                                    // 遍历 L 内部的所有 Entry (E_i)
                                    for (auto j = leafnode->children.begin(); j != leafnode->children.end(); ++j) {
                                        
                                        Entry* entry_child = dynamic_cast<Entry*>(j->second);
                                        if (!entry_child) continue;

										// #pragma omp critical (GlobalLoggerLock)
										// globalLogger.accessCandidate("candidate"); // 统计 *可能* 的 candidate

                                        // --- 5. 阶段一：廉价的 "Multilateration" 过滤 ---
                                        
                                        // d(Anchor, Candidate) (即 d(P_L, E_i))
                                        double d_anchor_candidate = entry_child->distanceToParent;

                                        // 应用三角不等式过滤： |d(Q,P_L) - d(E_i,P_L)| <= range
                                        if (std::abs(d_anchor_query - d_anchor_candidate) <= _query->range) {
                                            
                                            // --- 6. 阶段二：昂贵的 "精细过滤" (Final Check) ---
                                            // 只有通过了廉价过滤的 Entry 才计算真实距离
                                            
                                            // #pragma omp critical (GlobalLoggerLock)
                                            // globalLogger.accessCandidate("candidate"); // 统计 *可能* 的 candidate

                                            double final_distance = _query->_mtree->distance_function(_query->data, entry_child->data);

                                            if (final_distance <= _query->range) {
                                                // 确认是真实匹配

												#pragma omp critical (GlobalLoggerLock)
												globalLogger.accessCandidate("candidate"); // 统计 *可能* 的 candidate
                                                nearestQueue.push({entry_child, final_distance, final_distance});
                                            }
                                        }
                                    } // 结束 Entry 循环
                                } else {
                                    // 这是一个 InternalNode，将其推入队列继续遍历
                                    Node* internalNode = dynamic_cast<Node*>(child);
                                    assert(internalNode != NULL);
                                    pendingQueue.push({internalNode, childDistance, childMinDistance});
                                }
                            }
                        }
                    } // 结束 children 循环

                    if(pendingQueue.empty()) {
                        nextPendingMinDistance = std::numeric_limits<double>::infinity();
                    } else {
                        nextPendingMinDistance = pendingQueue.top().minDistance;
                    }

					// for(typename Node::ChildrenMap::const_iterator i = node->children.begin(); i != node->children.end(); ++i) {
					// 	IndexItem* child = i->second;
					// 	// --- 1. 第一次剪枝检查 (基于父节点到子节点距离) ---
					// 	if(std::abs(pending.distance - child->distanceToParent) - child->radius <= _query->range) {
					// 		double childDistance = _query->_mtree->distance_function(_query->data, child->data);
					// 		double childMinDistance = std::max(childDistance - child->radius, 0.0);
					// 		// --- 2. 第二次剪枝检查 (基于最小距离) ---
					// 		if(childMinDistance <= _query->range) {
					// 			LeafNode* leafnode = dynamic_cast<LeafNode*>(child);
					// 			// if (leafnode != NULL) {
					// 			// 	_query->nodeAccess++; // 计数 LeafNode 的 node access
					// 			// 	_query->accessedLeafNodeRadii.push_back(leafnode->radius);

					// 			// 	// 遍历该 LeafNode 的所有子 Entry，将它们全部推入 nearestQueue
					// 			// 	for (typename LeafNode::ChildrenMap::const_iterator j = leafnode->children.begin(); j != leafnode->children.end(); ++j) {
					// 			// 		// LeafNode 的子节点都是 Entry。
					// 			// 		Entry* entry_child = (Entry*)j->second; // 假设 Entry* 获取方式
					// 			// 		// 这里的距离使用 LeafNode 的父节点距离或不计算，因为目标是返回所有 Entry
					// 			// 		// 为了保持 result_item 的结构完整，我们暂时使用叶节点中心距离作为 Entry 距离
					// 			// 		double entryDistance = childDistance; 
										
					// 			// 		#pragma omp critical (GlobalLoggerLock)
					// 			// 		globalLogger.accessCandidate("candidate"); // 计数 Entry

					// 			// 		// 将 Entry 推入 nearestQueue。这里不需要进行距离检查
					// 			// 		nearestQueue.push({entry_child, entryDistance, 0.0});
					// 			// 	}
					// 			// } else {
					// 			// 	// 如果不是叶子节点，它必须是一个内部 Node
					// 			// 	Node* internalNode = dynamic_cast<Node*>(child);
					// 			// 	assert(internalNode != NULL);
									
					// 			// 	// 将内部节点推入 pendingQueue，以便稍后展开
					// 			// 	pendingQueue.push({internalNode, childDistance, childMinDistance});
					// 			// }
					// 			if(leafnode != NULL) {
					// 				//叶子节点
					// 				_query->nodeAccess++; // 计数叶子节点的 node access
					// 				_query->accessedLeafNodeRadii.push_back(leafnode->radius);
					// 			}
					// 			Entry* entry = dynamic_cast<Entry*>(child);
					// 			if(entry != NULL) { // 数据节点
					// 				#pragma omp critical (GlobalLoggerLock)
					// 				globalLogger.accessCandidate("candidate"); // 计数 Entry
					// 				nearestQueue.push({entry, childDistance, childMinDistance});
					// 			} else {
					// 				Node* node = dynamic_cast<Node*>(child);
					// 				assert(node != NULL);
					// 				pendingQueue.push({node, childDistance, childMinDistance});
					// 			}
					// 		}
					// 	}
					// }

					// if(pendingQueue.empty()) {
					// 	nextPendingMinDistance = std::numeric_limits<double>::infinity();
					// } else {
					// 	nextPendingMinDistance = pendingQueue.top().minDistance;
					// }

				}

				isEnd = true;
			}


			bool prepareNextNearest() {
				if(!nearestQueue.empty()) {
					ItemWithDistances<Entry> nextNearest = nearestQueue.top();
					if(nextNearest.distance <= nextPendingMinDistance) {
						nearestQueue.pop();
						currentResultItem.data = nextNearest.item->data;
						currentResultItem.distance = nextNearest.distance;
						++yieldedCount;
						return true;
					}
				}

				return false;
			}


			const query* _query;
			result_item currentResultItem;
			bool isEnd;
			std::priority_queue<ItemWithDistances<Node>> pendingQueue;
			double nextPendingMinDistance;
			std::priority_queue<ItemWithDistances<Entry>> nearestQueue;
			size_t yieldedCount;
		};


		/**
		 * @brief Begins the execution of the query and returns an interator
		 *        which refers to the first result.
		 */
		iterator begin() const {
			return iterator(this);
		}


		/**
		 * @brief Returns an iterator which informs that there are no more
		 *        results.
		 */
		iterator end() const {
			return {};
		}

	private:
		const mtree* _mtree;
		Data data;
		double range;
		size_t limit;
		friend class iterator;
	};



	enum {
		/**
		 * @brief The default minimum capacity of nodes in an M-Tree, when not
		 * specified in the constructor call.
		 */
		DEFAULT_MIN_NODE_CAPACITY = 50
	};

public:
	struct range_query_result {
        /** @brief 在范围内找到的项的列表。 */
        std::vector<typename query::result_item> matches;
        
        /** @brief 搜索期间访问的节点数。 */
        size_t nodeAccesses;

		/** @brief 搜索期间访问的 LeafNode 的半径列表。 */
        std::vector<double> leafNodeRadii;
    };


	/**
	 * @brief The main constructor of an M-Tree.
	 *
	 * @param min_node_capacity The minimum capacity of the nodes of an M-Tree.
	 *        Should be at least 2.
	 * @param max_node_capacity The maximum capacity of the nodes of an M-Tree.
	 *        Should be greater than @c min_node_capacity. If -1 is passed, then
	 *        the value <code>2*min_node_capacity - 1</code> is used.
	 * @param distance_function An instance of @c DistanceFunction.
	 * @param split_function An instance of @c SplitFunction.
	 *
	 * @param leaf_radius_threshold The threshold of leaf node, when radius < threshold
	 * 		  the min_node_capacity can be overlooked.
	 */
	explicit mtree(
			size_t min_node_capacity = DEFAULT_MIN_NODE_CAPACITY,
			size_t max_node_capacity = -1,
			double leaf_radius_threshold = 5, //叶子节点半径允许的阈值
			size_t min_compact_thresh = 5, // 触发紧凑性分裂的最小节点容量
			double pre_split_ratio = 2,     // 触发紧凑性分裂的半径膨胀比例
			const DistanceFunction& distance_function = DistanceFunction(),
			const SplitStrategy& split_strategy = SplitStrategy()
			// const SplitFunction& split_function = SplitFunction()
		)
		: minNodeCapacity(min_node_capacity),
		  maxNodeCapacity(max_node_capacity),
		  leafRadiusThreshold(leaf_radius_threshold),
		  minCompactnessThreshold(min_compact_thresh),
		  preSplitRadiusRatio(pre_split_ratio),
		  root(NULL),
		  distance_function(distance_function),
		  split_strategy(split_strategy)
		//   split_function(split_function)
	{
		if(max_node_capacity == size_t(-1)) {
			this->maxNodeCapacity = 2 * min_node_capacity - 1;
		}
	}

	// Cannot copy!
	mtree(const mtree&) = delete;

	// ... but moving is ok.
	/** @brief Move constructor. */
	mtree(mtree&& that)
		: root(that.root),
		  maxNodeCapacity(that.maxNodeCapacity),
		  minNodeCapacity(that.minNodeCapacity),
		  distance_function(that.distance_function),
		  split_strategy(that.split_strategy)
		//   split_function(that.split_function)
	{
		that.root = NULL;
	}


	~mtree() {
		delete root;
	}

	// Cannot copy!
	mtree& operator=(const mtree&) = delete;

	// ... but moving is ok.
	/** @brief Move assignment. */
	mtree& operator=(mtree&& that) {
		if(&that != this) {
			std::swap(this->root, that.root);
			this->minNodeCapacity = that.minNodeCapacity;
			this->maxNodeCapacity = that.maxNodeCapacity;
			this->distance_function = std::move(that.distance_function);
			this->split_strategy = std::move(that.split_strategy);
			// this->split_function = std::move(that.split_function);
		}
		return *this;
	}


	/**
	 * @brief Adds and indexes a data object.
	 * @details An object that is already indexed should not be added. There is
	 *          no validation, and the behavior is undefined if done.
	 * @param data The data object to index.
	 */
	void add(const Data& data) {
		if(root == NULL) {
			root = new RootLeafNode(data);
			root->addData(data, 0, this);
		} else {
			double distance = distance_function(data, root->data);
			try {
				root->addData(data, distance, this);
			} catch(SplitNodeReplacement& e) {
				Node* newRoot = new RootNode(root->data);
				delete root;
				root = newRoot;
				for(int i = 0; i < e.newNodes.size(); ++i) {
					Node* newNode = e.newNodes[i];
					double distance = distance_function(root->data, newNode->data);
					root->addChild(newNode, distance, this);
				}
			} catch(CompactnessPromotion& promo) {
				std::cout << "root node compact!" << std::endl;
				// --- CATCH 2: 紧凑性提升 (新增) ---
				// 根节点拒绝了 promo.promotedNode
				// 我们需要创建一个新的根节点来同时容纳旧根节点和这个新节点
				
				Node* newRoot = new RootNode(Data()); // 新的根
				Node* oldRoot = root;                     // 旧的根
				Node* promotedNode = promo.promotedNode; // 被提升的节点
				
				root = newRoot; // 更新 mtree 的根

				// 1. 将旧根节点添加为新根的子节点
				double oldRootDist = distance_function(root->data, oldRoot->data);
				root->addChild(oldRoot, oldRootDist, this);

				// 2. 将被提升的节点添加为新根的子节点
				double promotedDist = distance_function(root->data, promotedNode->data);
				root->addChild(promotedNode, promotedDist, this);
			}
		}
	}


	/**
	 * @brief Removes a data object from the M-Tree.
	 * @param data The data object to be removed.
	 * @return @c true if and only if the object was found.
	 */
	bool remove(const Data& data) {
		if(root == NULL) {
			return false;
		}

		double distanceToRoot = distance_function(data, root->data);
		try {
			root->removeData(data, distanceToRoot, this);
		} catch(RootNodeReplacement& e) {
			delete root;
			root = e.newRoot;
		} catch(DataNotFound) {
			return false;
		}
		return true;
	}


	/**
	 * @brief Performs a nearest-neighbors query on the M-Tree, constrained by
	 *        distance.
	 * @param query_data The query data object.
	 * @param range The maximum distance from @c query_data to fetched neighbors.
	 * @return A @c query object.
	 */
	range_query_result get_nearest_by_range(const Data& query_data, double range) const {
		// 1. 创建 lazy_query 对象
        query lazy_query = get_nearest(query_data, range, std::numeric_limits<unsigned int>::max());
        
        // 2. 准备结果结构体
        range_query_result result;
        result.nodeAccesses = 0;
        
        // 3. 通过迭代执行查询。
        //    填充 lazy_query.nodeAccesses
        for (const auto& item : lazy_query) {
            result.matches.push_back(item);
        }
        
        // 4. 从 query 对象中检索统计数据
        result.nodeAccesses = lazy_query.nodeAccess;
		result.leafNodeRadii = std::move(lazy_query.accessedLeafNodeRadii);

		// log
		globalLogger.accessIndex("mtree node");
        
        return result;
	}

	/**
	 * @brief Performs a nearest-neighbors query on the M-Tree, constrained by
	 *        the number of neighbors.
	 * @param query_data The query data object.
	 * @param limit The maximum number of neighbors to fetch.
	 * @return A @c query object.
	 */
	query get_nearest_by_limit(const Data& query_data, size_t limit) const {
		return get_nearest(query_data, std::numeric_limits<double>::infinity(), limit);
	}

	/**
	 * @brief Performs a nearest-neighbor query on the M-Tree, constrained by
	 *        distance and/or the number of neighbors.
	 * @param query_data The query data object.
	 * @param range The maximum distance from @c query_data to fetched neighbors.
	 * @param limit The maximum number of neighbors to fetch.
	 * @return A @c query object.
	 */
	query get_nearest(const Data& query_data, double range, size_t limit) const {
		return {this, query_data, range, limit};
	}

	/**
	 * @brief Performs a nearest-neighbor query on the M-Tree, without
	 *        constraints.
	 * @param query_data The query data object.
	 * @return A @c query object.
	 */
	query get_nearest(const Data& query_data) const {
		return {
			this,
			query_data,
			std::numeric_limits<double>::infinity(),
			std::numeric_limits<unsigned int>::max()
		};
	}

    size_t size() const {
        return root->size();
    }

public:

	void _check() const {
#ifndef NDEBUG
			if(root != NULL) {
				root->_check(this);
			}
#endif
	}

public:

using LayerNodeCounts = std::map<size_t, size_t>; // <层级, 节点数量>
LayerNodeCounts get_structure_info() const {
    // 结构信息：BFS 遍历
    if (root == NULL) {
        return {};
    }

    LayerNodeCounts counts;
    // 使用 std::queue 存储要遍历的节点指针
    std::queue<Node*> q; 
    
    // 根节点是第 0 层（或者第 1 层，取决于您的定义，这里用 0 方便数组索引）
    size_t current_layer = 0;

    // 存储当前层节点的数量（作为分隔符）
    q.push(root);
    q.push(nullptr); // 使用 nullptr 作为层级分隔符

    size_t node_count = 0;

    while (!q.empty()) {
        Node* node = q.front();
        q.pop();

        if (node == nullptr) {
            // 遇到层级分隔符
            if (node_count > 0) {
                // 记录上一层的节点数
                counts[current_layer] = node_count; 
                
                // 推进到下一层
                current_layer++;
                node_count = 0;
                
                // 如果队列中还有元素（即下一层有节点），则放入下一个分隔符
                if (!q.empty()) {
                    q.push(nullptr); 
                }
            }
            continue;
        }

        // 统计当前层节点
        node_count++; 

        // 将当前节点的所有子节点（下一层）加入队列
        for (auto const& [key, val] : node->children) {
            Node* child_node = dynamic_cast<Node*>(val);
            // 只有非叶节点（InternalNode 或 RootNode）才有 Node* 类型的子节点
            if (child_node != NULL) { 
                q.push(child_node);
            }
        }
    }
    
    // 最后的叶节点层（Layer N）的计数在循环中被记录，无需额外处理。
    return counts;
}

// LayerOverlapCounts 存储每一层（父层）对下一层（子层）的重叠统计
struct LayerOverlapInfo {
    size_t total_child_nodes;        // 下一层 (h+1) 的总枢轴数
    size_t multi_contained_children; // 下一层中被 > 1 个父节点覆盖的枢轴数
    size_t total_overlap_count;      // 下一层所有枢轴被覆盖的总次数 (Sigma(count_i) - total_child_nodes)
};

using LayerOverlapResults = std::map<size_t, LayerOverlapInfo>; // <层级 h, 统计信息>
// 在 mtree 类定义中添加此方法

LayerOverlapResults get_overlap_info() const {
    if (root == NULL) {
        return {};
    }

    LayerOverlapResults overlap_results;
    std::queue<Node*> q; 
    q.push(root);
    
    size_t current_layer = 0; // 当前层是父节点层 (N_i)

    while (!q.empty()) {
        size_t level_size = q.size();
        std::vector<Node*> current_level_nodes;

        // 1. 提取当前层的所有父节点 N_i
        for (size_t i = 0; i < level_size; ++i) {
            current_level_nodes.push_back(q.front());
            q.pop();
        }
        
        // 2. 收集下一层 (h+1) 所有子节点/数据项 C 的枢轴 P_C
        // next_level_containment_counts: <枢轴指针, 被包含次数>
        std::map<const Data*, size_t, typename functions::DataPtrLess<Data>> next_level_containment_counts;
        std::vector<Node*> next_level_nodes; // 用于 BFS 推进

        for (Node* parent_node : current_level_nodes) {
            for (auto const& [key, val] : parent_node->children) {
                const IndexItem* child_item = val; 
                const Data* pivot_ptr = &(child_item->data);

                // 确保每个枢轴只被初始化一次
                if (next_level_containment_counts.find(pivot_ptr) == next_level_containment_counts.end()) {
                     next_level_containment_counts[pivot_ptr] = 0;
                }
                
                Node* child_node = dynamic_cast<Node*>(val);
                if (child_node != NULL) {
                    next_level_nodes.push_back(child_node);
                }
            }
        }
        
        size_t total_next_level_pivots = next_level_containment_counts.size();
        if (total_next_level_pivots == 0) break;

        // 3. 计算子节点被同层父节点包含的次数 (O(N_h * N_{h+1}))
        for (Node* parent_node : current_level_nodes) {
            
            const double parent_radius = parent_node->radius;
            const Data& parent_pivot = parent_node->data;

            for (auto& pair : next_level_containment_counts) {
                const Data& child_pivot = *(pair.first);
                size_t& count = pair.second; // 引用，用于更新

                double distance = distance_function(parent_pivot, child_pivot);

                // 如果 N_i 的覆盖球包含子节点 C 的枢轴 P_C
                if (distance <= parent_radius) {
                    count++; 
                }
            }
        }
        
        // 4. 统计重叠数量和总次数，并记录结果
        size_t multi_contained_children = 0; // 被 > 1 次包含的枢轴数量
        size_t total_overlap_count = 0;      // 总重叠次数 (所有多余的包含)

        for (const auto& pair : next_level_containment_counts) {
            // pair.second 是该枢轴被同层节点包含的总次数 (count)
            if (pair.second > 1) { 
                multi_contained_children++;
                // 贡献的重叠次数 = (总次数 - 1)
                total_overlap_count += (pair.second - 1); 
            }
        }

        // 记录结果
        overlap_results[current_layer] = {
            total_next_level_pivots,
            multi_contained_children,
            total_overlap_count
        };

        // 5. 推进到下一层
        if (!next_level_nodes.empty()) {
            for (Node* next_node : next_level_nodes) {
                q.push(next_node);
            }
            current_layer++;
        } else {
            break;
        }
    }

    return overlap_results;
}

// 在 mtree 公有部分添加

void print_overlap_info() const {
    LayerOverlapResults results = get_overlap_info();

    if (results.empty()) {
        std::cout << "M-Tree 重叠度信息为空 (树结构未达到多层)." << std::endl;
        return;
    }

    std::cout << "--- M-Tree 重叠度分析 (子节点多路径包含) ---" << std::endl;
    for (const auto& pair : results) {
        size_t layer_index = pair.first;
        const auto& info = pair.second;
        
        double overlap_ratio = (info.total_child_nodes > 0) 
                               ? (double)info.multi_contained_children / info.total_child_nodes
                               : 0.0;
        
        std::cout << "父层级: " << layer_index 
                  << " | 子节点总数: " << info.total_child_nodes
                  << " | 被多次覆盖枢轴数: " << info.multi_contained_children
                  << " (" << std::fixed << std::setprecision(2) << (overlap_ratio * 100.0) << "%)"
                  << " | 总重叠次数: " << info.total_overlap_count
                  << std::endl;
    }
}
	// // typedef std::pair<Data, Data> PromotedPair;
	// typedef std::set<Data> Partition;
	// typedef std::vector<std::pair<Data, Partition>> SplitResult;


	size_t minNodeCapacity;
	size_t maxNodeCapacity;
	size_t leafRadiusThreshold;
	Node* root;

protected:
	DistanceFunction distance_function;
	SplitStrategy split_strategy;
	// SplitFunction split_function;

public:
	class IndexItem {
	public:
		Data data;
		double radius;
		double distanceToParent;

		virtual ~IndexItem() { };

		IndexItem() = delete;
		IndexItem(const IndexItem&) = delete;
		IndexItem(IndexItem&&) = delete;
		IndexItem& operator=(const IndexItem&) = delete;
		IndexItem& operator=(IndexItem&&) = delete;

	protected:
		IndexItem(const Data& data)
			: data(data),
			  radius(0),
			  distanceToParent(-1)
			{ }

	public:
		virtual size_t _check(const mtree* mtree) const {
			_checkRadius();
			_checkDistanceToParent();
			return 1;
		}

        virtual size_t size() const {
	    //printf("I am here...%zu\n", sizeof(Data) + 2 * sizeof(double));
            return sizeof(Data) + 2 * sizeof(double);
        }

	private:
		void _checkRadius() const {
			assert(radius >= 0);
		}

	protected:
		virtual void _checkDistanceToParent() const {
			assert(dynamic_cast<const RootNodeTrait*>(this) == NULL);
			assert(distanceToParent >= 0);
		}
	};


private:
	class Node : public IndexItem {
	public:
		virtual ~Node() {
			for(typename ChildrenMap::iterator i = children.begin(); i != children.end(); ++i) {
				IndexItem* child = i->second;
				delete child;
			}
		}

		virtual Node* newNode(const mtree* mtree) const = 0;

		void addData(const Data& data, double distance, const mtree* mtree) {
			doAddData(data, distance, mtree);
			checkMaxCapacity(mtree);
		}

        virtual size_t size() const override {
            auto fanout = DEFAULT_MIN_NODE_CAPACITY * 2 - 1;
            //auto fanout = children.size();
            auto res = sizeof(IndexItem) + fanout * (sizeof(Data) + sizeof(IndexItem*));
            for (auto& child: children) {
                res += child.second->size();
            }
	    //printf("I am a node: %zu %zu\n", children.size(), res);
            return res;
        }

#ifndef NDEBUG
		size_t _check(const mtree* mtree) const override {
			IndexItem::_check(mtree);
			_checkMinCapacity(mtree);
			_checkMaxCapacity(mtree);

			bool   childHeightKnown = false;
			size_t childHeight = 0;
			for(typename ChildrenMap::const_iterator i = children.begin(); i != children.end(); ++i) {
#ifndef NDEBUG
				const Data& data = i->first;
#endif
				IndexItem* child = i->second;

				assert(child->data == data);
				_checkChildClass(child);
				_checkChildMetrics(child, mtree);

				size_t height = child->_check(mtree);
				if(childHeightKnown) {
					assert(childHeight == height);
				} else {
					childHeight = height;
					childHeightKnown = true;
				}
			}

			return childHeight + 1;
		}
#endif

		typedef std::map<Data, IndexItem*> ChildrenMap;

		ChildrenMap children;

	protected:
		Node(const Data& data) : IndexItem(data) { }

		Node() : IndexItem(*((Data*)(0))) { assert(!"THIS SHOULD NEVER BE CALLED"); };

		Node(const Node&) = delete;
		Node(Node&&) = delete;
		Node& operator=(const Node&) = delete;
		Node& operator=(Node&&) = delete;

		virtual void doAddData(const Data& data, double distance, const mtree* mtree) = 0;

		virtual void doRemoveData(const Data& data, double distance, const mtree* mtree) = 0;

	public:
		void checkMaxCapacity(const mtree* mtree) {
			if(children.size() > mtree->maxNodeCapacity) {

				// 判断是否是 LeafNode
				const LeafNode* leafNode = dynamic_cast<const LeafNode*>(this);
				// 只有当它是 LeafNode 时，才应用半径超载规则
			double currentRadius = this->radius;
       		if (leafNode != nullptr) {
					if (currentRadius < mtree->leafRadiusThreshold) {
						// 尽管超载，但是半径紧凑，可以跳过分裂
						return;
					}
				}
				Partition all_entries;
				for(typename ChildrenMap::iterator i = children.begin(); i != children.end(); ++i) {
					all_entries.insert(i->first);
				}

				cached_distance_function_type cachedDistanceFunction(mtree->distance_function);

				// 更新后的 partition 方法
				SplitResult<Data> partitions_result = mtree->split_strategy(all_entries, cachedDistanceFunction, currentRadius);
				// 检查分裂结果是否有效（至少要有两个分区）
				assert(partitions_result.size() >= 2);

				// 创建一个用于替代当前节点的新节点列表，长度为 k
				std::vector<Node*> newNodes;
				
				// 遍历 SplitResult 中的每个分区
				for (const auto& partition_pair : partitions_result) {
					const Data& promotedData = partition_pair.first;
					const Partition& partition = partition_pair.second;

					// a. 为当前分区创建新的节点 (Node)
					// newSplitNodeReplacement(promotedData) 应该返回一个新的 InternalNode 或 LeafNode
					Node* newNode = newSplitNodeReplacement(promotedData);
					
					// b. 将原节点的孩子（IndexItem*）移动到新节点中
					for(typename Partition::iterator j = partition.begin(); j != partition.end(); ++j) {
						const Data& data = *j; // 这是孩子节点的 key (Entry key)
						
						// 找到原节点中对应的孩子指针
						// 注意：由于 all_entries 已被 splitStrategy 清空，我们需要使用 children 映射
						// 确保孩子节点在 children 中存在
						assert(children.count(data)); 
						
						IndexItem* child = children[data]; // 获取孩子指针
						children.erase(data);              // 从当前节点移除孩子

						// 计算新推广数据到孩子的距离
						double distance = cachedDistanceFunction(promotedData, data);
						
						// 将孩子添加到新的节点中
						newNode->addChild(child, distance, mtree);
					}
					
					// 将新创建的节点添加到列表中
					newNodes.push_back(newNode);
				}

				// Partition secondPartition;
				// PromotedPair promoted = mtree->split_function(all_entries, secondPartition, cachedDistanceFunction);

				// Node* newNodes[2];
				// for(int i = 0; i < 2; ++i) {
				// 	Data& promotedData    = (i == 0) ? promoted.first : promoted.second;
				// 	Partition& partition = (i == 0) ? all_entries : secondPartition;

				// 	Node* newNode = newSplitNodeReplacement(promotedData);
				// 	for(typename Partition::iterator j = partition.begin(); j != partition.end(); ++j) {
				// 		const Data& data = *j;
				// 		IndexItem* child = children[data];
				// 		children.erase(data);
				// 		double distance = cachedDistanceFunction(promotedData, data);
				// 		newNode->addChild(child, distance, mtree);
				// 	}

				// 	newNodes[i] = newNode;
				// }
				assert(children.empty());

				throw SplitNodeReplacement(newNodes);
			}

		}

	protected:
		virtual Node* newSplitNodeReplacement(const Data&) const = 0;

	public:
		virtual void addChild(IndexItem* child, double distance, const mtree* mtree) = 0;

		virtual void removeData(const Data& data, double distance, const mtree* mtree) {
			doRemoveData(data, distance, mtree);
			if(children.size() < getMinCapacity(mtree)) {
				throw NodeUnderCapacity();
			}
		}

		virtual size_t getMinCapacity(const mtree* mtree) const = 0;

	public:
		void updateMetrics(IndexItem* child, double distance) {
			child->distanceToParent = distance;
			this->updateRadius(child);
		}

		void updateRadius(IndexItem* child) {
			this->radius = std::max(this->radius, child->distanceToParent + child->radius);
		}


		virtual void _checkMinCapacity(const mtree* mtree) const = 0;

	private:
		void _checkMaxCapacity(const mtree* mtree) const {
			assert(children.size() <= mtree->maxNodeCapacity);
		}

	protected:
		virtual void _checkChildClass(IndexItem* child) const = 0;

	private:
#ifndef NDEBUG
		void _checkChildMetrics(IndexItem* child, const mtree* mtree) const {
			double dist = mtree->distance_function(child->data, this->data);
			assert(child->distanceToParent == dist);

			/* TODO: investigate why the following line
			 * 		assert(child->distanceToParent + child->radius <= this->radius);
			 * is not the same as the code below:
			 */
			double sum = child->distanceToParent + child->radius;
			assert(sum <= this->radius);
		}
#endif
	};


	class RootNodeTrait : public virtual Node {
		void _checkDistanceToParent() const {
			assert(this->distanceToParent == -1);
		}
	};


	class NonRootNodeTrait : public virtual Node {
		size_t getMinCapacity(const mtree* mtree) const override {
			return mtree->minNodeCapacity;
		}

		void _checkMinCapacity(const mtree* mtree) const override {
			assert(this->children.size() >= mtree->minNodeCapacity);// 
		}
	};


	class LeafNodeTrait : public virtual Node {
		void doAddData(const Data& data, double distance, const mtree* mtree) {
			// 1. 紧凑性检查 (在叶节点级别)
			// (distance 是新数据点到当前叶节点中心的距离)
			if (this->children.size() >= mtree->minCompactnessThreshold &&
				distance > (this->radius * mtree->preSplitRadiusRatio)) 
			{
				// std::cout << "leaf node compact!" << std::endl;
				// --- 失败：触发紧凑性提升 ---
				// "新建一个节点，抛出异常
				Node* promotedNode = this->newNode(mtree);
				promotedNode->data = data;
				// 使用 addChild 插入，避免递归检查
				promotedNode->addChild(new Entry(data), 0.0, mtree); 
				
				// 抛出异常，阻止数据插入到当前节点
				throw CompactnessPromotion(promotedNode);
			}
			
			// 2. 紧凑性检查通过：正常插入
			Entry* entry = new Entry(data);
			assert(this->children.find(data) == this->children.end());
			this->children[data] = entry;
			assert(this->children.find(data) != this->children.end());
			this->updateMetrics(entry, distance);
		}

		void addChild(IndexItem* child, double distance, const mtree* mtree) {
			assert(this->children.find(child->data) == this->children.end());
			this->children[child->data] = child;
			assert(this->children.find(child->data) != this->children.end());
			this->updateMetrics(child, distance);
		}

		Node* newSplitNodeReplacement(const Data& data) const {
			return new LeafNode(data);
		}

		void doRemoveData(const Data& data, double distance, const mtree* mtree) {
			if(this->children.erase(data) == 0) {
				throw DataNotFound{data};
			}
		}


		void _checkChildClass(IndexItem* child) const {
			assert(dynamic_cast<Entry*>(child) != NULL);
		}
	};


	class NonLeafNodeTrait : public virtual Node {
		void doAddData(const Data& data, double distance, const mtree* mtree) {
			struct CandidateChild {
				Node* node;
				double distance;
				double metric;
			};

			CandidateChild minRadiusIncreaseNeeded = { NULL, -1.0, std::numeric_limits<double>::infinity() };
			CandidateChild nearestDistance         = { NULL, -1.0, std::numeric_limits<double>::infinity() };

			for(typename Node::ChildrenMap::iterator i = this->children.begin(); i != this->children.end(); ++i) {
				Node* child = dynamic_cast<Node*>(i->second);
				assert(child != NULL);
				double distance = mtree->distance_function(child->data, data);
				if(distance > child->radius) {
					double radiusIncrease = distance - child->radius;
					if(radiusIncrease < minRadiusIncreaseNeeded.metric) {
						minRadiusIncreaseNeeded = { child, distance, radiusIncrease };
					}
				} else {
					if(distance < nearestDistance.metric) {
						nearestDistance = { child, distance, distance };
					}
				}
			}

			CandidateChild chosen = (nearestDistance.node != NULL)
			                      ? nearestDistance
			                      : minRadiusIncreaseNeeded;

			Node* child = chosen.node;
			try {
				child->addData(data, chosen.distance, mtree);
				this->updateRadius(child);
			} catch(CompactnessPromotion& promo) {
				// std::cout << "non leaf node compact!" << std::endl;
				// --- 场景 1: 捕获到紧凑性提升 (C1 抛出了 C_n) ---
				Node* promotedChildNode = promo.promotedNode; // 这就是 C_n
				// "判断 B1到 C_n的距离会不会导致B1 的半径扩大太大"
				double newDistance = mtree->distance_function(this->data, promotedChildNode->data);

				// 3. 父节点 (B1) 检查自己的紧凑性
				if (this->children.size() >= mtree->minCompactnessThreshold &&
					newDistance > (this->radius * mtree->preSplitRadiusRatio))
				{
					// --- 失败：B1 也被破坏了 ---
					// "新建一个B_n，其子节点是C_n，然后继续向上抛出异常"
					// 创建 B_n
					Node* newParentNode = this->newNode(mtree);
					// 将 C_n 作为 B_n 的子节点
					newParentNode->addChild(promotedChildNode, newDistance, mtree);
					// 将 B_n 向上抛给 A
					throw CompactnessPromotion(newParentNode);
				} else {
					this->addChild(promotedChildNode, newDistance, mtree);
				}
			} catch(SplitNodeReplacement& e) {
				// Replace current child with new nodes
#ifndef NDEBUG
				size_t _ =
#endif
					this->children.erase(child->data);
				assert(_ == 1);
				delete child;

				for(int i = 0; i < e.newNodes.size(); ++i) {
					Node* newChild = e.newNodes[i];
					double distance = mtree->distance_function(this->data, newChild->data);
					addChild(newChild, distance, mtree);
				}
			}
		}


		void addChild(IndexItem* newChild_, double distance, const mtree* mtree) {
			Node* newChild = dynamic_cast<Node*>(newChild_);
			assert(newChild != NULL);

			struct ChildWithDistance {
				Node* child;
				double distance;
			};

			std::vector<ChildWithDistance> newChildren;
			newChildren.push_back(ChildWithDistance{newChild, distance});

			while(!newChildren.empty()) {
				ChildWithDistance cwd = newChildren.back();
				newChildren.pop_back();

				newChild = cwd.child;
				distance = cwd.distance;
				typename Node::ChildrenMap::iterator i = this->children.find(newChild->data);
				if(i == this->children.end()) {
					this->children[newChild->data] = newChild;
					this->updateMetrics(newChild, distance);
				} else {
					Node* existingChild = dynamic_cast<Node*>(this->children[newChild->data]);
					assert(existingChild != NULL);
					assert(existingChild->data == newChild->data);

					// Transfer the _children_ of the newChild to the existingChild
					for(typename Node::ChildrenMap::iterator i = newChild->children.begin(); i != newChild->children.end(); ++i) {
						IndexItem* grandchild = i->second;
						double distanceToExistingChild = mtree->distance_function(existingChild->data, grandchild->data);
						existingChild->addChild(grandchild, distanceToExistingChild, mtree);
					}
					newChild->children.clear();
					delete newChild;

					this->updateRadius(existingChild);

					try {
						existingChild->checkMaxCapacity(mtree);
					} catch(SplitNodeReplacement& e) {
#ifndef NDEBUG
						size_t _ =
#endif
							this->children.erase(existingChild->data);
						assert(_ == 1);
						delete existingChild;

						for(int i = 0; i < e.newNodes.size(); ++i) {
							Node* newNode = e.newNodes[i];
							double distance = mtree->distance_function(this->data, newNode->data);
							newChildren.push_back(ChildWithDistance{newNode, distance});
						}
					}
				}
			}
		}


		Node* newSplitNodeReplacement(const Data& data) const {
			return new InternalNode(data);
		}


		void doRemoveData(const Data& data, double distance, const mtree* mtree) {
			for(typename Node::ChildrenMap::iterator i = this->children.begin(); i != this->children.end(); ++i) {
				Node* child = dynamic_cast<Node*>(i->second);
				assert(child != NULL);
				if(std::abs(distance - child->distanceToParent) <= child->radius) {
					double distanceToChild = mtree->distance_function(data, child->data);
					if(distanceToChild <= child->radius) {
						try {
							child->removeData(data, distanceToChild, mtree);
							this->updateRadius(child);
							return;
						} catch(DataNotFound&) {
							// If DataNotFound was thrown, then the data was not found in the child
						} catch(NodeUnderCapacity&) {
							Node* expandedChild = balanceChildren(child, mtree);
							this->updateRadius(expandedChild);
							return;
						}
					}
				}
			}

			throw DataNotFound{data};
		}


		Node* balanceChildren(Node* theChild, const mtree* mtree) {
			// Tries to find anotherChild which can donate a grand-child to theChild.

			Node* nearestDonor = NULL;
			double distanceNearestDonor = std::numeric_limits<double>::infinity();

			Node* nearestMergeCandidate = NULL;
			double distanceNearestMergeCandidate = std::numeric_limits<double>::infinity();

			for(typename Node::ChildrenMap::iterator i = this->children.begin(); i != this->children.end(); ++i) {
				Node* anotherChild = dynamic_cast<Node*>(i->second);
				assert(anotherChild != NULL);
				if(anotherChild == theChild) continue;

				double distance = mtree->distance_function(theChild->data, anotherChild->data);
				if(anotherChild->children.size() > anotherChild->getMinCapacity(mtree)) {
					if(distance < distanceNearestDonor) {
						distanceNearestDonor = distance;
						nearestDonor = anotherChild;
					}
				} else {
					if(distance < distanceNearestMergeCandidate) {
						distanceNearestMergeCandidate = distance;
						nearestMergeCandidate = anotherChild;
					}
				}
			}

			if(nearestDonor == NULL) {
				// Merge
				for(typename Node::ChildrenMap::iterator i = theChild->children.begin(); i != theChild->children.end(); ++i) {
					IndexItem* grandchild = i->second;
					double distance = mtree->distance_function(grandchild->data, nearestMergeCandidate->data);
					nearestMergeCandidate->addChild(grandchild, distance, mtree);
				}

				theChild->children.clear();
				this->children.erase(theChild->data);
				delete theChild;
				return nearestMergeCandidate;
			} else {
				// Donate
				// Look for the nearest grandchild
				IndexItem* nearestGrandchild = nullptr;
				double nearestGrandchildDistance = std::numeric_limits<double>::infinity();
				for(typename Node::ChildrenMap::iterator i = nearestDonor->children.begin(); i != nearestDonor->children.end(); ++i) {
					IndexItem* grandchild = i->second;
					double distance = mtree->distance_function(grandchild->data, theChild->data);
					if(distance < nearestGrandchildDistance) {
						nearestGrandchildDistance = distance;
						nearestGrandchild = grandchild;
					}
				}

#ifndef NDEBUG
				size_t _ =
#endif
					nearestDonor->children.erase(nearestGrandchild->data);
				assert(_ == 1);
				theChild->addChild(nearestGrandchild, nearestGrandchildDistance, mtree);
				return theChild;
			}
		}


		void _checkChildClass(IndexItem* child) const {
			assert(dynamic_cast<InternalNode*>(child) != NULL
			   ||  dynamic_cast<LeafNode*>(child)     != NULL);
		}
	};


	class RootLeafNode : public RootNodeTrait, public LeafNodeTrait {
	public:
		RootLeafNode(const Data& data) : Node(data) { }

		void removeData(const Data& data, double distance, const mtree* mtree) override {
			try {
				Node::removeData(data, distance, mtree);
			} catch (NodeUnderCapacity&) {
				assert(this->children.empty());
				throw RootNodeReplacement{NULL};
			}
		}

		size_t getMinCapacity(const mtree* mtree) const override {
			return 1;
		}

		void _checkMinCapacity(const mtree* mtree) const override {
			assert(this->children.size() >= 1);
		}

		RootLeafNode(const mtree* mtree) : Node(Data()) {}

		Node* newNode(const mtree* mtree) const override {
			return new RootLeafNode(mtree);
		}
	};

	class RootNode : public RootNodeTrait, public NonLeafNodeTrait {
	public:
		RootNode(const Data& data) : Node(data) {}

		RootNode(const mtree* mtree) : Node(Data()) {}
		
		Node* newNode(const mtree* mtree) const override {
			return new RootNode(mtree);
		}

	private:
		void removeData(const Data& data, double distance, const mtree* mtree) override {
			try {
				Node::removeData(data, distance, mtree);
			} catch(NodeUnderCapacity&) {
				// Promote the only child to root
				Node* theChild = dynamic_cast<Node*>(this->children.begin()->second);
				Node* newRoot;
				if(dynamic_cast<InternalNode*>(theChild) != NULL) {
					newRoot = new RootNode(theChild->data);
				} else {
					assert(dynamic_cast<LeafNode*>(theChild) != NULL);
					newRoot = new RootLeafNode(theChild->data);
				}

				for(typename Node::ChildrenMap::iterator i = theChild->children.begin(); i != theChild->children.end(); ++i) {
					IndexItem* grandchild = i->second;
					double distance = mtree->distance_function(newRoot->data, grandchild->data);
					newRoot->addChild(grandchild, distance, mtree);
				}
				theChild->children.clear();

				throw RootNodeReplacement{newRoot};
			}
		}


		size_t getMinCapacity(const mtree* mtree) const override {
			return 2;
		}

		void _checkMinCapacity(const mtree* mtree) const override {
			assert(this->children.size() >= 2);
		}
	};


	class InternalNode : public NonRootNodeTrait, public NonLeafNodeTrait {
	public:
		InternalNode(const Data& data) : Node(data) { }

		InternalNode(const mtree* mtree) : Node(Data()) {}

		Node* newNode(const mtree* mtree) const override {
			return new InternalNode(mtree);
		}
	};


	class LeafNode : public NonRootNodeTrait, public LeafNodeTrait {
	public:
		LeafNode(const Data& data) : Node(data) { }

		LeafNode(const mtree* mtree) : Node(Data()) {}

		Node* newNode(const mtree* mtree) const override {
			return new LeafNode(mtree);
		}
	};


	class Entry : public IndexItem {
	public:
		Entry(const Data& data) : IndexItem(data) { }

        virtual size_t size() const {
            return 0;
        }
	};
};



} /* namespace mtree */



#endif /* MTREE_H_ */
