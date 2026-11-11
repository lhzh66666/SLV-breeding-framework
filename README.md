# DGSR Calculation and Visualization Code

This repository contains the R code and processed/raw data used for the **Molecular Horticulture** manuscript. It includes the modular validation of DGSR components (Figure 2), the multi-environment phenotypic data (Figure 3), and the cross-species end-to-end proof-of-concept evaluation using the soybean RIL population (Figure 4).

## File Structure and Description:

### Figure 2 (Modular Validation)
*   `Fig2A_data.csv`: Processed data for Panel A (A0 vs. Proportion of stable QTLs).
*   `Fig2B_data.csv`: Processed data for Panel B (G0 vs. Observed F1 progeny performance).
*   `Fig2C_data.csv`: Processed data for Panel C (P0 vs. Environmental sensitivity).
*   `README_Figure2.md`: Detailed instructions for the Figure 2 dataset.

### Figure 3 (Phenotypic Data)
*   `Data.xlsx`: Raw phenotypic data across multiple environments.
*   `A_ge_results.xlsx`: Additive × Environment (A×E) interaction effects.

### Figure 4 (End-to-End Soybean Validation)
*   `01_DGSR_Calculation.R`: Script to calculate G0 (selection index), A0 (favorable alleles), P0 (polymorphism), and the composite DGSR score.
*   `02_Figure_4_Visualization.R`: Script to generate the 3-panel composite plot (Figure 4A, 4B, 4C).
*   `Soybean_DGSR_EndToEnd.csv`: The output file containing DGSR scores for the 118 RILs.
*   `MAPMAKER.txt`: Genotypic data (259 polymorphic markers) for the soybean RIL population.

## Requirements:
R version 4.0 or higher. The required packages are `readxl`, `ggplot2`, `ggrepel`, `patchwork`, `dplyr`, and `scales`.

## Usage:
For Figure 4 reproducibility, ensure `Data.xlsx` and `MAPMAKER.txt` are in the same working directory.
1. Run `01_DGSR_Calculation.R` to generate the DGSR scores.
2. Run `02_Figure_4_Visualization.R` to generate the high-resolution Figure 4.