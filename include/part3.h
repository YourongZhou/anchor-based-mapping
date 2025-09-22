// include/part3.h
#pragma once
#include "anchor_common.h"

extern "C" {
    // inputs:
    //   anchor_map_path : 文件，格式 anchor_id \t bound_distance\n
    //   reference_fasta : 要检索的 reference 库
    //   out_candidates_path : 输出候选 ID（每行一个 ID），若为空则 stdout
    //
    // 语义（矬矬版）:
    //   逐条 reference r，计算 dist(r, anchor_i) 对所有 anchor_i；
    //   如果对每个 anchor_i 都有 dist(r,anchor_i) <= bound_i ，则把 r 作为 candidate。
    int retrieve_candidates(const char* anchor_map_path, const char* reference_fasta, const char* out_candidates_path);
}
