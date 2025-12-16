# Anchor-Based Mapping
A high-performance C++ implementation of three approximate sequence matching strategies for genomic data: M-Tree, Seed-and-Extend, and Anchor-based indexing methods.

## Features
Three Retrieval Strategies: M-Tree metric space indexing, Seed-and-Extend with FM-Index, and Anchor-based pre-computed distance tables
Comprehensive Evaluation: Built-in metrics calculation (precision, recall, access counts) and performance benchmarking
Parallel Processing: OpenMP-based parallelization for query processing
Flexible Configuration: Extensive command-line parameters for method tuning
Testing Framework: Python-based testing infrastructure with Jupyter notebooks for visualization
## Quick Start
Build
```bash
# Build the main executable 
make anchor_mapping
```
Basic Usage
```bash
# M-Tree method (default)  
./anchor_mapping --method mtree --maxDist 3 --num_queries 100  
  
# Seed-and-Extend method  
./anchor_mapping --method sae --maxDist 3 --seed_len 16 --num_queries 100  
  
# Anchor-based method  
./anchor_mapping --method anchor --maxDist 3 --anchor_len 250 --num_anchors 1000 --num_queries 100
```
Command Line Options
| Parameter      | Description                                              | Default |
|----------------|----------------------------------------------------------|---------|
| --method       | Retrieval strategy (mtree, sae, anchor)                  | mtree   |
| --maxDist      | Maximum edit distance (comma-separated list supported)   | 3       |
| --anchor_len   | Length of anchor sequences                               | 150      |
| --num_anchors  | Number of anchors for anchor method                     | 100     |
| --seed_len     | Seed length for SAE method                               | 10      |

## Architecture
The system implements three distinct retrieval strategies with a unified evaluation framework:

### M-Tree Strategy
Hierarchical metric tree indexing with triangle inequality pruning
Extracts k-mers from reference and builds searchable tree structure
Returns positions, node access counts, and leaf node radii tools.h:61-64
### Seed-and-Extend Strategy
FM-Index construction for compressed reference storage
Exact seed matching with gapped X-drop extension
Dynamic seed length calculation based on query parameters tools.h:67-73
### Anchor-Based Strategy
Pre-computed anchor distance tables
Hash-based lookups with set intersection
Tunable anchor selection and distance filtering candidate_retriever.cpp:219-223

## Testing and Evaluation
The project includes a comprehensive testing framework using Jupyter notebooks:

```bash
cd test/  
jupyter notebook figures_upd.ipynb  # Main evaluation notebook
```

### Key Metrics
- Position-based: True Positives, False Positives, False Negatives
- Distance-based: Average and maximum edit distances
- Access patterns: Index and candidate access counts
- Performance: Precision, Recall, FP/TP ratios figures_upd.ipynb:37-48
### Batch Experiment Execution
The testing framework supports batch parameter sweeps:

```python
# Example from testing notebook  
maxDist_values = [1, 3, 5, 7, 9, 11]  
df_sae = run_batch_experiment('sae', maxDist_values, fixed_seed_len='16')  
df_mtree = run_batch_experiment('mtree', maxDist_values)  
df_anchor = run_batch_experiment('anchor', maxDist_values)
figures_upd.ipynb:121-134
```

## Dependencies
- SeqAn2: Sequence analysis library for FM-Index and sequence operations
- OpenMP: Parallel processing support
- Python 3: Testing framework (pandas, matplotlib, seaborn)
- C++17: Modern C++ features
  
## File Structure
> 
├── main.cpp                 # Main executable entry point  
├── include/  
│   ├── tools.h             # Core function declarations  
│   ├── metrics.h           # Metrics calculation functions  
│   └── data_types.h        # Data structure definitions  
├── src/  
│   ├── tools.cpp           # Core implementations  
│   └── candidate_retriever.cpp  # Anchor-based retrieval  
├── test/  
│   ├── figures_upd.ipynb   # Main evaluation notebook  
│   ├── figures_upd_sae.ipynb  # SAE-specific tests  
│   └── figures_upd_anchor.ipynb  # Anchor-specific tests  
└── index_cache/            # Cached index files  
## Example Output
===== Average per-query Metrics =====  
Average TP: 1  
Average FP: 11.54  
Average FN: 0  
Average Recall: 1  
Average Precision: 0.0797448  
Average FP/TP: 11.54  
Average average distance: 10.1049  
Average maximum distance: 14.26  
=== Access Summary ===  
Index accesses:     600  
Candidate accesses: 1665  
Total accesses:     2265  
output_sae.txt:129-142

## Notes
The system supports batch execution of multiple parameter values in a single run for efficient testing
Ground truth generation uses parallel processing for scalability
All three methods share a common metrics evaluation framework for fair comparison
The testing framework generates publication-quality visualizations including precision/recall curves, access count charts, and distance distribution plots main.cpp:352-412
I deeply appreciate your interest.:)

