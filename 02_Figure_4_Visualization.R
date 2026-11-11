# =====================================================================
# Script 2: Figure 4 Visualization
# Description: Generates the 3-panel composite plot (A, B, C).
# Input: Soybean_DGSR_EndToEnd.csv
# Output: Figure_4_Final.pdf/png
# =====================================================================

# 1. Load Packages
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(scales)

# 2. Load Data
df <- read.csv("Soybean_DGSR_EndToEnd.csv")
df_clean <- df %>% filter(!is.na(DGSR)) 
top10 <- df_clean %>% arrange(desc(DGSR)) %>% slice(1:10)
df_clean$Is_Top10 <- ifelse(df_clean$Line_ID %in% top10$Line_ID, "Top 10", "Others")

# 3. Panel A: Bubble Plot (A0 vs G0)
pA <- ggplot(df_clean, aes(x = A0, y = G0_proxy)) +
  geom_point(aes(size = P0, color = DGSR), alpha = 0.7) +
  geom_point(data = top10, shape = 24, size = 3, fill = "red", color = "black") +
  geom_text_repel(data = top10, aes(label = Line_ID), size = 3, max.overlaps = 20) +
  scale_color_gradientn(colors = c("#F0F8FF", "#4682B4", "#FF4500"), name = "DGSR") +
  scale_size_continuous(range = c(1, 5), name = "P0") +
  theme_classic() +
  labs(x = "A0 (Genetic Potential)", y = "G0 (Phenotypic Performance)", 
       title = "A  DGSR Evaluation of 114 Soybean RILs") +
  theme(plot.title = element_text(face = "bold", size = 12))

# 4. Panel B: Waterfall Plot (DGSR Ranking)
pB <- ggplot(df_clean, aes(x = reorder(Line_ID, -DGSR), y = DGSR, fill = Is_Top10)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("Others" = "grey70", "Top 10" = "#E64B35"), name = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "114 Soybean RILs (Ranked by DGSR)", y = "DGSR Value", 
       title = "B  End-to-End Ranking of DGSR")

# 5. Panel C: Correlation Plot (DGSR vs Yield)
cor_res <- cor.test(df_clean$DGSR, df_clean$YLD, method = "pearson")
r_val <- round(cor_res$estimate, 3)
p_val <- format.pval(cor_res$p.value, digits = 3)

pC <- ggplot(df_clean, aes(x = DGSR, y = YLD)) +
  geom_point(aes(color = Is_Top10), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed") +
  geom_text_repel(data = top10, aes(label = Line_ID), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Others" = "#4DBBD5", "Top 10" = "#E64B35"), name = "") +
  theme_classic() +
  annotate("text", x = min(df_clean$DGSR)*1.1, y = max(df_clean$YLD)*0.9, 
           label = paste0("R = ", r_val, "\np = ", p_val), hjust = 0, size = 4) +
  labs(x = "DGSR Value", y = "Yield (YLD)", 
       title = "C  Correlation between DGSR and Yield")

# 6. Combine and Save
final_plot <- pA / (pB | pC) + plot_layout(heights = c(1, 1))
ggsave("Figure_4_Final.pdf", final_plot, width = 12, height = 8, dpi = 300)
ggsave("Figure_4_Final.png", final_plot, width = 12, height = 8, dpi = 300)
print("Figure 4 generated successfully!")