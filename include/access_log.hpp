#ifndef ACCESS_LOGGER_HPP
#define ACCESS_LOGGER_HPP

#include <iostream>
#include <string>
#include <atomic>
#include <mutex>

class AccessLogger {
public:
    // 构造函数
    AccessLogger() : index_access(0), candidate_access(0) {}

    // 记录一次 index 访问
    inline void accessIndex(const std::string& key) {
        {
            std::lock_guard<std::mutex> lock(mu);
            ++index_access;
        }
#ifdef ACCESSLOGGER_VERBOSE
        std::cout << "[Index Access] key=" << key << std::endl;
#endif
    }

    // 记录一次 candidate 访问
    inline void accessCandidate(const std::string& cid) {
        {
            std::lock_guard<std::mutex> lock(mu);
            ++candidate_access;
        }
#ifdef ACCESSLOGGER_VERBOSE
        std::cout << "[Candidate Access] cid=" << cid << std::endl;
#endif
    }

    // 输出统计结果
    void report() const {
        std::cout << "=== Access Summary ===" << std::endl;
        std::cout << "Index accesses:     " << index_access.load() << std::endl;
        std::cout << "Candidate accesses: " << candidate_access.load() << std::endl;
        std::cout << "Total accesses:     " 
                  << (index_access.load() + candidate_access.load()) 
                  << std::endl;
    }

    // 清零
    void reset() {
        index_access.store(0);
        candidate_access.store(0);
    }

private:
    std::atomic<size_t> index_access;
    std::atomic<size_t> candidate_access;
    mutable std::mutex mu;
};

// 声明全局访问日志对象
extern AccessLogger globalLogger;

#endif // ACCESS_LOGGER_HPP