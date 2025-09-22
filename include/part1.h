// include/part1.h
#pragma once
#include "anchor_common.h"
#include <string>

extern "C" {
    // 从 input_fasta 随机抽取 num_anchors 条序列，写入 out_anchor_fasta (FASTA)
    // 返回 0 成功，非 0 失败
    int generate_anchors(const char* input_fasta, const char* out_anchor_fasta, int num_anchors, unsigned int seed);
}
