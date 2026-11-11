# =====================================================================
# Script 1: End-to-end DGSR Calculation Pipeline (Soybean RIL Panel)
# Description: Calculates G0, A0, P0, and composite DGSR score.
# Input: Data.xlsx (phenotype), MAPMAKER.txt (genotype)
# Output: Soybean_DGSR_EndToEnd.csv
# =====================================================================

# 1. Setup and Load Packages
library(readxl)
setwd("your_working_directory_path") # Please change to your local path

# 2. Load Phenotypic Data
pheno_raw <- read_excel("Data.xlsx", skip = 1)
colnames(pheno_raw) <- c("LL","LW","LSL","LC","PN","SWPP1","PH","ST","HP","SNN","BN",
                         "SPPN1","SPPN2","SPPN3","SPPN4","LDG","FD","ETM","PL2","PB2",
                         "PCT","SWPP2","WP","YLD","SL","SB","SW","SWPP3","SCN","PRO","OIL")
pheno <- pheno_raw[1:118, ]
pheno$Line_ID <- paste0("RIL", sprintf("%03d", 1:118))

# 3. Calculate G0 (Selection Index: YLD 0.5, PRO 0.25, OIL 0.25)
pheno_scaled <- scale(pheno[, c("YLD", "PRO", "OIL")])
g0_raw <- pheno_scaled[, "YLD"] * 0.5 + pheno_scaled[, "PRO"] * 0.25 + pheno_scaled[, "OIL"] * 0.25
pheno$G0_proxy <- (g0_raw - min(g0_raw, na.rm=TRUE)) / (max(g0_raw, na.rm=TRUE) - min(g0_raw, na.rm=TRUE))

# 4. Load Genotypic Data (MAPMAKER format)
genotype_lines <- readLines("MAPMAKER.txt")
marker_lines <- grep("^\\*", genotype_lines, value = TRUE)
marker_lines <- marker_lines[1:259]  # Keep 259 polymorphic markers
geno_matrix <- do.call(rbind, strsplit(marker_lines, "\\s+"))
geno_data <- geno_matrix[, -1] # Keep 118 families
geno_num <- matrix(ifelse(geno_data == "A", 1, ifelse(geno_data == "B", 0, NA)), nrow = nrow(geno_data))

# 5. Calculate A0 (Accumulation of favorable alleles)
A0_score <- rep(0, 118) 
for (i in 1:nrow(geno_num)) { # Iterate over markers (rows)
  geno_i <- geno_num[i, ] 
  valid_idx <- !is.na(geno_i) & !is.na(pheno$YLD) 
  
  if(sum(valid_idx) > 50) {
    y_A <- pheno$YLD[valid_idx & geno_i == 1]
    y_B <- pheno$YLD[valid_idx & geno_i == 0]
    
    if(length(y_A) > 0 && length(y_B) > 0) {
      mean_A <- mean(y_A, na.rm = TRUE)
      mean_B <- mean(y_B, na.rm = TRUE)
      
      if(is.finite(mean_A) && is.finite(mean_B) && mean_A > mean_B) {
        A0_score[valid_idx & geno_i == 1] <- A0_score[valid_idx & geno_i == 1] + 1
      }
    }
  }
}
pheno$A0 <- A0_score / max(A0_score, na.rm = TRUE)

# 6. Calculate P0 (Polymorphism and representativeness)
P0_score <- colSums(geno_num == 1, na.rm = TRUE) * 0.5 + colSums(geno_num == 0, na.rm = TRUE) * 0.5
pheno$P0 <- (P0_score - min(P0_score)) / (max(P0_score) - min(P0_score))

# 7. Calculate final DGSR
pheno$DGSR <- (pheno$A0 + ifelse(is.na(pheno$G0_proxy), 0, pheno$G0_proxy) + pheno$P0) / 10 * 0.1

# 8. Save Results
write.csv(pheno[, c("Line_ID", "A0", "G0_proxy", "P0", "DGSR", "YLD", "PRO", "OIL")], 
          "Soybean_DGSR_EndToEnd.csv", row.names = FALSE)
print("DGSR Calculation Finished!")