# ==========================================================
# Opsonophagocytosis Data Clenaing and Analysis
# ==========================================================

# 1) Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)

# 2) Source processing functions
source("opsono_data_function.R")

# 3) Load raw data
df <- read.csv("raw_data.csv")

# 4) Clean and structure data
df_clean <- parse_opsono_data(df)

# 5) Compute Salmonella clearance
df_clearance <- compute_salmonella_clearance(df_clean)

# 6) Run ANOVA + Pairwise Post Hoc Tests
res_clearance <- run_genotype_anova(df_clearance, "clearance_pct")
print(res_clearance$anova_table)
print(res_clearance$pairwise_matrix)

res_mrc1_ssc <- run_genotype_anova(df_clean, "mrc1_salpos_mean_ssc")
print(res_mrc1_ssc$anova_table)
print(res_mrc1_ssc$pairwise_matrix)

res_salpos_pct <- run_genotype_anova(df_clean, "salpos_pct")
print(res_salpos_pct$anova_table)
print(res_salpos_pct$pairwise_matrix)

res_salpos_fitc <- run_genotype_anova(df_clean, "salpos_mfi_fitc")
print(res_salpos_fitc$anova_table)
print(res_salpos_fitc$pairwise_matrix)

res_salpos_pe <- run_genotype_anova(df_clean, "salpos_mfi_pe")
print(res_salpos_pe$anova_table)
print(res_salpos_pe$pairwise_matrix)

res_salneg_pe <- run_genotype_anova(df_clean, "salneg_mfi_pe")
print(res_salneg_pe$anova_table)
print(res_salneg_pe$pairwise_matrix)

res_salpos_amcyan <- run_genotype_anova(df_clean, "salpos_mfi_amcyan")
print(res_salpos_amcyan$anova_table)
print(res_salpos_amcyan$pairwise_matrix)

res_salneg_amcyan <- run_genotype_anova(df_clean, "salneg_mfi_amcyan")
print(res_salneg_amcyan$anova_table)
print(res_salneg_amcyan$pairwise_matrix)

res_total_pe <- run_genotype_anova(df_clean, "total_mfi_pe")
print(res_total_pe$anova_table)
print(res_total_pe$pairwise_matrix)

res_total_amcyan <- run_genotype_anova(df_clean, "total_mfi_amcyan")
print(res_total_amcyan$anova_table)
print(res_total_amcyan$pairwise_matrix)

# 6) Visualize data
# Salmonell clearance
plot_genotype_bargraph(df_clearance, "clearance_pct")
# MFI / Uptake Comparisons
plot_genotype_bargraph(df_clean, "mrc1_salpos_mean_ssc")
plot_genotype_bargraph(df_clean, "salpos_pct")
plot_genotype_bargraph(df_clean, "salpos_mfi_fitc")
plot_genotype_bargraph(df_clean, "salpos_mfi_pe")
plot_genotype_bargraph(df_clean, "salneg_mfi_pe")
plot_genotype_bargraph(df_clean, "salpos_mfi_amcyan")
plot_genotype_bargraph(df_clean, "salneg_mfi_amcyan")
plot_genotype_bargraph(df_clean, "total_mfi_pe")
plot_genotype_bargraph(df_clean, "total_mfi_amcyan")