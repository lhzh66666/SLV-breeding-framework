====================================================================
DGSR (Degree of Germplasm Self-Renewal) Framework
====================================================================
This repository contains the R scripts and processed data files used to 
generate Figure 2 (Modular validation of DGSR components in soybean) 
for the manuscript: "Degree of Germplasm Self-Renewal (DGSR): A Multi-Omics-Informed 
Quantitative Framework for Plant Genetic Resource Phenotyping and Breeding Prioritization".

The modular validation independently tests the biological interpretability 
of the three DGSR components: superior allele inventory (A0), realized 
genetic contribution (G0), and diversity-exchange potential (P0).

--------------------------------------------------------------------
1. SYSTEM REQUIREMENTS
--------------------------------------------------------------------
- Software: R (version 4.0.0 or higher) and RStudio
- Required R Packages:
  ggplot2
  ggpubr
  patchwork

To install the required packages in R, run:
install.packages(c("ggplot2", "ggpubr", "patchwork"))

--------------------------------------------------------------------
2. FILE STRUCTURE
--------------------------------------------------------------------
./data/
  - Fig2A_data.csv : Processed data for Panel A (A0 vs. Proportion of stable QTLs).
  - Fig2B_data.csv : Processed data for Panel B (G0 vs. Observed F1 progeny performance).
  - Fig2C_data.csv : Processed data for Panel C (P0 vs. Environmental sensitivity).

./scripts/
  - Figure2_plotting_script.R : R script to perform Spearman correlation analysis 
                                and generate the multi-panel figure (Figure 2).

./outputs/
  - Figure_2_Complete.tiff    : Final generated figure (300 dpi).
  - Figure_2_Complete.pdf     : Vector version of the final figure.

--------------------------------------------------------------------
3. HOW TO RUN (REPRODUCIBILITY)
--------------------------------------------------------------------
Step 1: Download this repository to your local machine.
Step 2: Open RStudio and set the working directory to the root folder 
        of this repository (e.g., setwd("path/to/DGSR-framework")).
Step 3: Ensure the CSV files are located in the ./data/ folder.
Step 4: Open and run the script:
        source("scripts/Figure2_plotting_script.R")
Step 5: The script will calculate the Spearman rank correlation coefficients (rho) 
        and save the generated figures into the ./outputs/ folder.

--------------------------------------------------------------------
4. DATA AVAILABILITY
--------------------------------------------------------------------
The data provided here are the processed summary statistics extracted from 
the independent soybean validation datasets (118 QTLs, diallel-cross data, 
and ecophysiological records). The full raw phenotypic and genotypic datasets 
are provided in Supplementary Tables S1 and S2 of the manuscript.

--------------------------------------------------------------------
5. CONTACT
--------------------------------------------------------------------
For any questions, bug reports, or data inquiries, please contact:
Huizhen Liang
Institute of Chinese Herbal Medicines, Henan Academy of Agricultural Sciences
Email: lhzh66666@163.com
ORCID: https://orcid.org/0000-0002-9709-0721