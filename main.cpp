// main.cpp
#include "experiment_manager.h"
#include <iostream>

int main(int argc, char* argv[]) {
    ExperimentConfig config;

    // ----- 解析命令行参数 -----
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--truncate_ref_len" && i + 1 < argc) {
            config.truncate_ref_len = std::stoul(argv[++i]);
        } else if (arg == "--anchor_len" && i + 1 < argc) {
            config.anchor_len = std::stoul(argv[++i]);
        } else if (arg == "--num_anchors" && i + 1 < argc) {
            config.num_anchors = std::stoi(argv[++i]);
        } else if (arg == "--num_queries" && i + 1 < argc) {
            config.num_queries = std::stoi(argv[++i]);
        } else if (arg == "--maxDist" && i + 1 < argc) {
            config.maxDist_list = parseIntList(argv[++i]);
        } else if (arg == "--use_all" && i + 1 < argc) {
            std::string val = argv[++i];
            config.use_all = (val == "true" || val == "1");
        } else if (arg == "--start_idx" && i + 1 < argc) {
            config.start_idx = std::stod(argv[++i]);
        } else if (arg == "--end_idx" && i + 1 < argc) {
            config.end_idx = std::stod(argv[++i]);
        } else if (arg == "--last_random" && i + 1 < argc) {
            std::string val = argv[++i];
            config.last_random = (val == "true" || val == "1");
        } else if (arg == "--seed" && i + 1 < argc) {
            config.seed = std::stoul(argv[++i]);
        } else if (arg == "--anchor_radius" && i + 1 < argc) {
            config.anchor_radius = std::stoi(argv[++i]);
        } else if (arg == "--min_node_capacity" && i + 1 < argc) {
            config.min_node_capacity = std::stoi(argv[++i]);
        } else if (arg == "--max_node_capacity" && i + 1 < argc) {
            config.max_node_capacity = std::stoi(argv[++i]);
        } else if (arg == "--leaf_radius_threshold" && i + 1 < argc) {
            config.leaf_radius_threshold = std::stod(argv[++i]);
        } else if (arg == "--compactness_min_capacity" && i + 1 < argc) {
            config.compactness_min_capacity = std::stoi(argv[++i]);
        } else if (arg == "--compactness_radius_factor" && i + 1 < argc) {
            config.compactness_radius_factor = std::stod(argv[++i]);
        } else if (arg == "--split_strategy" && i + 1 < argc) {
            config.split_strategy_name = argv[++i];
        } else if (arg == "--auto_gen" && i + 1 < argc) {
            std::string val = argv[++i];
            config.auto_gen = (val == "true" || val == "1");
        } else if (arg == "--seed_len" && i + 1 < argc) {
            config.seed_len = std::stoi(argv[++i]);
        } else if (arg == "--method" && i + 1 < argc) {
            config.method = argv[++i];
        } else if (arg == "--index_dir" && i + 1 < argc) {
            config.index_dir = argv[++i];
        } else if (arg == "--force_rebuild") {
            config.force_rebuild = true;
        } else if (arg == "--require_distance") {
            config.require_distance = true;
        }
    }

    try {
        ExperimentManager manager(config);
        manager.run();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
