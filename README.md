# SLV breeding framework analysis scripts

This repository contains the analysis scripts for the manuscript:

**A Systems Life View framework to quantify resource-allocation trade-offs and enable coordinated multi-trait improvement in specialized metabolite crops**

## Requirements
- Python 3.15
- pandas, numpy, scipy, matplotlib, openpyxl
- R 4.3.1 (for BLUP estimation, lme4 package)

## Files
- `SLV_analysis.py`: Main script for SLV model construction, ΔS calculation, time-split blind validation, LOYO cross-validation, Pareto front analysis, and permutation tests.
- `test_single.py`: Script for single-trait selection comparison.
- `Table1-4.xlsx`: Soybean phenotype and genotype data (Sheet 1: Phenotype; Sheet 2: Genotype).
- `Top_selection_results.xlsx`: Yield and isoflavone means, differences, and permutation test P values for Top10/20/30 selection intensities.
- `SLV_Detailed_Predictions.xlsx`: Predicted and observed ΔS for each RIL line.
- `A_ge_results.xlsx`: Genotype–environment adaptability A(g,e) for each line across locations.

## Usage
1. Ensure `Table1-4.xlsx` is in the same directory as `SLV_analysis.py`.
2. Run:
   ```python
   python SLV_analysis.py
   ```
3. Output files and figures will be generated in the current directory.

## Data
Soybean RIL genotype and phenotype data are also available in Supplementary Data S1 of the manuscript.

## Contact
Huizhen Liang (lhzh66666@163.com)