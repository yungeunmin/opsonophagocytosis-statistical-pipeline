# Opsonophagocytosis Statistical Pipeline

Modular R pipeline for analyzing flow cytometry–based opsonophagocytosis assays.

## Overview

This project provides a reproducible framework to:

- Parse FlowJo CSV exports  
- Extract experimental metadata from filenames  
- Compute Salmonella clearance  
- Run one-way ANOVA by genotype  
- Perform Tukey post-hoc comparisons  
- Generate bar plots with significance letters  

---

## Structure

opsonophagocytosis-statistical-pipeline/
│
├── sample_data/
│   └── raw_data_sample.csv  
│
├── R/
│   ├── opsono_data_function.R  
│   └── opsono_data_analysis.R  
│
└── README.md  

---

## Key Functions

- `parse_opsono_data()`  
  Cleans and structures raw FlowJo exports.

- `compute_salmonella_clearance()`  
  Calculates clearance:  
  (1 - (24hr / mean(2hr))) × 100  

- `run_genotype_anova()`  
  Performs ANOVA and Tukey HSD post-hoc testing.

- `plot_genotype_bargraph()`  
  Creates mean ± SE bar plots with compact letter display.

---

## Statistical Methods

- One-way ANOVA  
- Tukey HSD  
- Mean ± SE visualization  

---

## Reproducibility

All analyses are script-based.  
No manual spreadsheet steps required.