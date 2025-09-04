# Progresa-Oportunidades RDD Analysis - Data Preparation

## Overview

This repository contains the data preparation pipeline for analyzing Mexico's Progresa-Oportunidades conditional cash transfer program using regression discontinuity design (RDD). The analysis examines how the program's transition from Progresa (1997-2002) to Oportunidades (2002+) affected targeting and outcomes.

## Project Structure

```
issue-1.1-data-preparation/
│
├── code/                           # All R scripts
│   ├── 01_load_and_clean_data.R   # Load and clean raw data
│   ├── 02_merge_datasets.R        # Merge datasets
│   ├── 03_create_rdd_datasets.R   # Create RDD-specific datasets
│   └── helper_functions.R         # Reusable functions
│
├── data/
│   ├── raw/                       # Original data files (not tracked in git)
│   └── processed/                 # Clean datasets
│       ├── merged_entry_iml_rg.csv    # Main cross-sectional dataset
│       ├── iml_panel_90_10.csv        # Panel dataset
│       ├── regional_cutoffs.csv       # Optimal cutoffs by region
│       └── rdd_analysis_dataset.csv   # RDD-ready dataset
│
├── output/
│   ├── reports/                   # Data quality and diagnostic reports
│   │   ├── data_quality_report.txt
│   │   ├── merge_diagnostics.txt
│   │   └── data_preparation_summary.txt
│   │
│   └── tables/                    # Summary statistics
│       ├── summary_statistics.csv
│       └── balance_tests.csv
│
├── README.md                      # This file
└── CODEBOOK.md                    # Variable documentation
```

## Quick Start

### Prerequisites

Required R packages:
```r
install.packages(c(
  "tidyverse",    # Data manipulation
  "haven",        # Read Stata files
  "janitor",      # Data cleaning
  "rdrobust",     # RDD analysis
  "rddensity",    # Manipulation testing
  "here",         # Project paths
  "skimr",        # Data summaries
  "assertr",      # Validation
  "broom"         # Tidy outputs
))
```

### Running the Pipeline

1. **Set up your working directory**:
   ```r
   # Adjust this path to your local setup
   base_path <- "~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS"
   ```

2. **Run scripts in sequence**:
   ```r
   # Step 1: Load and clean raw data
   source("code/01_load_and_clean_data.R")
   
   # Step 2: Merge datasets
   source("code/02_merge_datasets.R")
   
   # Step 3: Create RDD datasets
   source("code/03_create_rdd_datasets.R")
   ```

3. **Or run all at once**:
   ```r
   # Source helper functions first
   source("code/helper_functions.R")
   
   # Run full pipeline
   source("code/01_load_and_clean_data.R")
   source("code/02_merge_datasets.R")
   source("code/03_create_rdd_datasets.R")
   ```

## Data Sources

### Required Raw Data Files

Place these files in your base directory:

1. **Marginalization Indices**:
   - `iml1995.dta` - 1995 marginalization index
   - `iml2000.dta` - 2000 marginalization index  
   - `iml2005.dta` - 2005 marginalization index

2. **Program Entry Data**:
   - `entry_loc_300710.dta` - Program entry years (source 1)
   - `entry_locality_enca.dta` - Program entry years (source 2)

3. **Regional Classification**:
   - `regiones41.dta` - Regional groupings of municipalities

4. **Panel Data**:
   - `Base_marginacion_localidad_90-10.csv` - Locality panel 1990-2010

## Key Outputs

### 1. Main Cross-sectional Dataset (`merged_entry_iml_rg.csv`)
- **Use for**: Main RDD analysis
- **Key variables**: `iml_95`, `post_2002`, `entry_year`, `region41`
- **N**: ~85,000 localities

### 2. Panel Dataset (`iml_panel_90_10.csv`)
- **Use for**: Temporal dynamics analysis
- **Key variables**: `mge_id`, `year`, `treated`, IML components
- **Structure**: Locality-year panel

### 3. Regional Cutoffs (`regional_cutoffs.csv`)
- **Use for**: Regional heterogeneity analysis
- **Key variables**: `region`, `optimal_cutoff`, `rd_estimate`
- **N**: 41 regions with valid cutoffs

### 4. RDD Analysis Dataset (`rdd_analysis_dataset.csv`)
- **Use for**: Normalized RDD analysis
- **Key variables**: `iml_95_normalized`, `above_cutoff`
- **Feature**: Running variable normalized by regional cutoffs

## Key Transformations

### 1. Locality Identification
- Created unique 9-digit locality ID: `EEMMMMLLLL`
- Standardized municipality codes to 3 digits

### 2. Program Entry
- Combined two entry data sources (priority: entry_30)
- Created binary indicators: `post_2002`, `post_2003`
- Categorized entry periods (Progresa vs Oportunidades)

### 3. RDD Variables
- Found optimal cutoffs by region using rdrobust
- Normalized running variable by regional cutoffs
- Created bandwidth-restricted samples

## Data Quality Checks

The pipeline includes extensive validation:

- Duplicate locality detection
- Missing value reporting
- Merge success rates
- Range validation for percentages
- Balance tests at cutoffs
- McCrary manipulation tests

## Next Steps

This data feeds into subsequent analyses:

- **Issue 1.2**: National RDD analysis
- **Issue 1.3**: Regional heterogeneity 
- **Issue 1.4**: Temporal dynamics
- **Issue 1.5**: Component analysis
- **Issue 1.6**: Robustness checks

## Troubleshooting

### Common Issues

1. **Encoding errors with panel data**:
   - Solution: File uses Latin-1 encoding, handled in script

2. **Memory issues with large merges**:
   - Solution: Scripts use efficient joins and filter early

3. **Missing regional cutoffs**:
   - Some regions have too few observations
   - Check `regional_cutoffs.csv` for valid regions

### Validation Checks

After running, verify:
```r
# Check output files exist
list.files("issue-1.1-data-preparation/data/processed/")

# Check dimensions
df <- read_csv("issue-1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")
dim(df)  # Should be ~85,000 rows

# Check treatment distribution
table(df$post_2002)
```

