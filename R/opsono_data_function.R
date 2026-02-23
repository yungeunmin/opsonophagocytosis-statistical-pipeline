library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(multcompView)

# -----------------------------------------------------------
# 1. Function: parse_opsono_data
#
# Purpose:
#   Take raw FlowJo CSV export and return a clean,
#   analysis-ready dataframe with:
#     - standardized measurement column names
#     - parsed experimental metadata from filename
#     - numeric conversion of percent columns
#
# Assumptions:
#   - Column structure matches current FlowJo export format
#   - Filename encodes genotype, serum concentration,
#     timepoint, treatment type, and replicate ID
# -----------------------------------------------------------

parse_opsono_data <- function(df_raw) {
  
  # -------------------------------------------------------
  # 1) Rename FlowJo measurement columns
  #
  # We replace long gating-path column names with
  # biologically meaningful variable names.
  #
  # Naming convention:
  #   salpos  = Salmonella-positive population
  #   salneg  = Salmonella-negative population
  #   mfi     = mean fluorescence intensity
  #   total   = total leucocyte parent population
  # -------------------------------------------------------
  
  df <- df_raw %>%
    rename(
      filename = X,  # FlowJo exports sample name as first unnamed column
      
      # SSC-A for mrc1+ within sal+ gate
      mrc1_salpos_mean_ssc = l_d.singles_f..singles_s.leu...Mean..SSC.A.,
      
      # Primary endpoint: percent Salmonella-positive cells
      salpos_pct = l_d.singles_f..singles_s.leu.sal....Freq..of.Parent,
      
      # MFI measurements in sal+ gate
      salpos_mfi_fitc = l_d.singles_f..singles_s.leu.sal....Mean..Comp.FITC.A.,
      salpos_mfi_pe = l_d.singles_f..singles_s.leu.sal....Mean..Comp.PE.A.,
      salpos_mfi_amcyan = l_d.singles_f..singles_s.leu.sal....Mean..Comp.AmCyan.A.,
      
      # MFI measurements in sal− gate
      salneg_mfi_pe = l_d.singles_f..singles_s.leu.sal....Mean..Comp.PE.A..1,
      salneg_mfi_amcyan = l_d.singles_f..singles_s.leu.sal....Mean..Comp.AmCyan.A..1,
      
      # Total leucocyte population MFIs
      total_mfi_pe = l_d.singles_f..singles_s.leu...Mean..Comp.PE.A.,
      total_mfi_amcyan = l_d.singles_f..singles_s.leu...Mean..Comp.AmCyan.A.
    )
  
  
  # -------------------------------------------------------
  # 2) Parse experimental design metadata from filename
  #
  # Example filename:
  #   opsonophagocytosis_TLR3_Serum(10ul)_2hr_seq_015.fcs
  #
  # Extract:
  #   genotype   → TLR3
  #   serum_ul   → 10
  #   timepoint  → 2hr
  #   treatment  → seq
  #   replicate  → 015
  # -------------------------------------------------------
  
  meta <- df %>%
    mutate(x_clean = str_remove(filename, "\\.fcs$")) %>%
    transmute(
      filename,
      
      # Knockout type or WT
      genotype  = str_extract(x_clean, "TLR[0-9]+|WT"),
      
      # Serum concentration in microliters
      serum_ul  = str_extract(x_clean, "(?<=Serum\\()\\d+(?=ul\\))") %>%
        as.numeric(),
      
      # Treatment time (2hr or 24hr)
      timepoint = str_extract(x_clean, "\\d+hr"),
      
      # Opsonization type (sequential or parallel)
      treatment = str_extract(x_clean, "seq|par"),
      
      # Replicate ID (final numeric block in filename)
      replicate = str_extract(x_clean, "(?<=_)\\d+$")
    )
  
  
  # -------------------------------------------------------
  # 3) Clean numeric measurement columns
  #
  # FlowJo exports percent columns as strings like:
  #   "12.4 %"
  #
  # We remove the percent sign and convert to numeric.
  #
  # All other MFI columns are already numeric.
  # -------------------------------------------------------
  
  df_clean <- df %>%
    mutate(
      salpos_pct = as.numeric(
        str_remove(salpos_pct, "\\s*%\\s*$")
      )
    ) %>%
    select(-filename)  # remove duplicate before binding
  
  
  # -------------------------------------------------------
  # 4) Combine metadata and measurement data
  #
  # Final dataframe contains:
  #   - experimental design variables
  #   - numeric measurement variables
  #   - ready for modeling or summary statistics
  # -------------------------------------------------------
  
  df_final <- bind_cols(meta, df_clean)
  
  return(df_final)
}

# -----------------------------------------------------------
# 2. Function: compute_salmonella_clearance
#
# Purpose:
#   Calculate Salmonella clearance (%) using thesis formula:
#
#   clearance = (1 - (%sal+ at 24hr / mean(%sal+ at 2hr))) * 100
#
# Grouping structure:
#   baseline mean at 2hr is computed within:
#     genotype × serum_ul × treatment
#
# Input:
#   df_clean (output of parse_opsono_data)
#
# Output:
#   Dataframe containing only 24hr rows with:
#     - baseline_2hr_mean
#     - clearance_pct
# -----------------------------------------------------------

compute_salmonella_clearance <- function(df_clean) {
  
  # -------------------------------
  # 1) Compute baseline (2hr mean)
  # -------------------------------
  
  baseline_2hr <- df_clean %>%
    filter(timepoint == "2hr") %>%
    group_by(genotype, serum_ul, treatment) %>%
    summarise(
      baseline_2hr_mean = mean(salpos_pct, na.rm = TRUE),
      n_2hr = sum(!is.na(salpos_pct)),
      .groups = "drop"
    )
  
  # -------------------------------
  # 2) Filter 24hr rows
  # -------------------------------
  
  df_24hr <- df_clean %>%
    filter(timepoint == "24hr")
  
  # -------------------------------
  # 3) Join baseline to 24hr rows
  # -------------------------------
  
  df_joined <- df_24hr %>%
    left_join(
      baseline_2hr,
      by = c("genotype", "serum_ul", "treatment")
    )
  
  # -------------------------------
  # 4) Compute clearance
  #
  # Handle edge cases:
  #   - baseline = 0 → NA
  #   - missing values → NA
  # -------------------------------
  
  df_clearance <- df_joined %>%
    mutate(
      clearance_pct = case_when(
        is.na(baseline_2hr_mean) ~ NA_real_,
        baseline_2hr_mean == 0 ~ NA_real_,
        is.na(salpos_pct) ~ NA_real_,
        TRUE ~ (1 - (salpos_pct / baseline_2hr_mean)) * 100
      )
    )
  
  return(df_clearance)
}

# -----------------------------------------------------------
# 3. Function: run_genotype_anova
#
# Purpose:
#   Run one-way ANOVA (outcome ~ genotype)
#   Perform Tukey post-hoc comparisons
#   Return ANOVA summary and pairwise p-value matrix
#
# Input:
#   data        - dataframe
#   outcome_var - string name of outcome column
#
# Output:
#   list containing:
#     - anova_table
#     - pairwise_matrix (symmetric matrix of adjusted p-values)
#     - tukey_raw
# -----------------------------------------------------------

run_genotype_anova <- function(data, outcome_var) {
  
  # ---------------------------------------
  # Ensure genotype is a factor
  # Order WT first if present
  # ---------------------------------------
  
  geno_levels <- unique(data$genotype)
  
  if ("WT" %in% geno_levels) {
    geno_levels <- c("WT", setdiff(geno_levels, "WT"))
  }
  
  data$genotype <- factor(data$genotype, levels = geno_levels)
  
  # ---------------------------------------
  # Build model formula dynamically
  # ---------------------------------------
  
  formula_obj <- as.formula(paste(outcome_var, "~ genotype"))
  
  # ---------------------------------------
  # Run ANOVA
  # ---------------------------------------
  
  model <- aov(formula_obj, data = data)
  anova_summary <- summary(model)
  
  # ---------------------------------------
  # Tukey post-hoc test
  # ---------------------------------------
  
  tukey_res <- TukeyHSD(model)$genotype
  
  comparisons <- rownames(tukey_res)
  pvals <- tukey_res[, "p adj"]
  
  # ---------------------------------------
  # Create symmetric p-value matrix
  # ---------------------------------------
  
  groups <- levels(data$genotype)
  
  pairwise_matrix <- matrix(
    NA,
    nrow = length(groups),
    ncol = length(groups),
    dimnames = list(groups, groups)
  )
  
  for (i in seq_along(comparisons)) {
    comps <- strsplit(comparisons[i], "-")[[1]]
    g1 <- comps[1]
    g2 <- comps[2]
    
    pairwise_matrix[g1, g2] <- pvals[i]
    pairwise_matrix[g2, g1] <- pvals[i]  # make symmetric
  }
  
  return(list(
    anova_table = anova_summary,
    pairwise_matrix = pairwise_matrix,
    tukey_raw = tukey_res
  ))
}


# -----------------------------------------------------------
# 4. Function: plot_genotype_bargraph
#
# Purpose:
#   Create barplot with:
#     - Mean ± SE
#     - ANOVA compact letter display
#
# Input:
#   data        - dataframe
#   outcome_var - string
#
# Output:
#   ggplot object
# -----------------------------------------------------------

plot_genotype_bargraph <- function(data, outcome_var) {
  
  # Ensure genotype is factor (WT first if present)
  geno_levels <- unique(data$genotype)
  if ("WT" %in% geno_levels) {
    geno_levels <- c("WT", setdiff(geno_levels, "WT"))
  }
  data$genotype <- factor(data$genotype, levels = geno_levels)
  
  # -------------------------------
  # Run ANOVA + Tukey
  # -------------------------------
  
  formula_obj <- as.formula(paste(outcome_var, "~ genotype"))
  model <- aov(formula_obj, data = data)
  
  tukey_res <- TukeyHSD(model)$genotype
  
  # Convert Tukey p-values to named vector
  tukey_p <- tukey_res[, "p adj"]
  names(tukey_p) <- rownames(tukey_res)
  
  # Generate compact letter display
  letters <- multcompView::multcompLetters(tukey_p)$Letters
  
  letters_df <- data.frame(
    genotype = names(letters),
    letters = letters
  )
  
  # -------------------------------
  # Summary statistics
  # -------------------------------
  
  summary_df <- data %>%
    group_by(genotype) %>%
    summarise(
      mean = mean(.data[[outcome_var]], na.rm = TRUE),
      se = sd(.data[[outcome_var]], na.rm = TRUE) /
        sqrt(sum(!is.na(.data[[outcome_var]]))),
      .groups = "drop"
    )
  
  summary_df <- left_join(summary_df, letters_df, by = "genotype")
  
  # -------------------------------
  # Plot
  # -------------------------------
  
  p <- ggplot(summary_df,
              aes(x = genotype, y = mean)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = mean - se,
                      ymax = mean + se),
                  width = 0.2) +
    geom_text(aes(label = letters,
                  y = mean + se + max(se, na.rm = TRUE) * 0.5),
              size = 5) +
    labs(
      x = "Genotype",
      y = outcome_var
    ) +
    theme_minimal()
  
  return(p)
}



