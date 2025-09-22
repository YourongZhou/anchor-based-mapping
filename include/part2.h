// include/part2.h
#pragma once
#include "anchor_common.h"

extern "C" {
    // 输入 query string，读取 anchors fasta，返回一个简单的文本文件 (anchor_id \t distance \n)
    // 如果 out_map_path 为空则打印到 stdout
    int compute_anchor_distances(const char* query_seq, const char* anchor_fasta, const char* out_map_path);
}
