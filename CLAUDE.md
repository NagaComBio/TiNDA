# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TiNDA (Tumor in Normal Detection Analysis) is an R package that rescues somatic variants misclassified as germline due to tumor DNA contamination in patient blood/control samples. It uses Canopy's EM-cluster algorithm to partition variants into clusters and classify them as somatic rescue, CHIP, or germline based on VAF patterns.

## Development Commands

```bash
# Run R CMD check
R CMD check .

# Install package locally
R CMD INSTALL .

# Build package
R CMD build .

# Run package tests
R -e "devtools::test()"

# Load package for development
devtools::load_all()

# Document with roxygen2
devtools::document()

# Run specific test file
R -e "testthat::test_file('tests/testthat/test-filename.R')"

# Check coverage
covr::package_coverage()
```

## Architecture

### Core Components

**R/TiNDA.R** - Main package logic
- `TiNDA()` function: Primary entry point that takes tumor/control variant data and classifies variants
- Uses Canopy's `canopy.cluster()` for EM-based clustering
- Key classification parameters: `max_control_af` (0.25), `min_tumor_af` (0.01), `min_clst_members` (0.85)
- Classification flow: Canopy clustering → potential somatic/CHIP cluster selection → germline exclusion → optional homozygous rescue

**R/data.R** - Data simulation utilities
- `generate_depth()`: Helper for generating variant depths with normal distribution
- `simulate_variants()`: Simulates germline, somatic, and CHIP variants with configurable VAFs
- `generate_test_data()`: Generates test data across all chromosomes

**R/plot.R** - Visualization functions
- `canopy_clst_plot()`: Shows Canopy cluster assignments in VAF space
- `tinda_clst_plot()`: Shows TiNDA final classification (Somatic_Rescue, Germline, CHIP)
- `tinda_linear_plot()`: Linear genome browser-style plot across chromosomes
- `tinda_summary_plot()`: Combined view with all plots and summary table

**R/write_data.R** - Output utilities
- `write_data()`: Exports TiNDA results to TSV

### Data Flow

1. Input: Data frame with CHR, POS, Control_ALT_DP, Control_DP, Tumor_ALT_DP, Tumor_DP
2. VAF calculation: Control_AF = Control_ALT_DP/Control_DP, Tumor_AF = Tumor_ALT_DP/Tumor_DP
3. Canopy EM clustering with 10 clusters
4. Cluster classification based on Area of Interest (AOI) criteria
5. Optional CHIP detection (controlled_af: 0.02-0.35, tumor_af < 0.25)
6. Homozygous TiN rescue at top-left corner
7. Output: Classified variants with TiN_Class column

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| max_control_af | 0.25 | Max control VAF for somatic rescue |
| min_tumor_af | 0.01 | Min tumor VAF for somatic rescue |
| min_clst_members | 0.85 | Min cluster members (proportion) in AOI |
| num_run | 1 | EM runs per cluster count |
| find_chip | TRUE | Enable CHIP detection |
| max_control_af_chip | 0.40 | Max control VAF for CHIP |
| max_tumor_af_chip | 0.25 | Max tumor VAF for CHIP |

### Dependencies

- Canopy (>= 1.3.0): EM clustering
- ggplot2: Visualization
- dplyr, purrr: Data manipulation
- readr: TSV I/O

### Supported Reference Genomes

- hg19: data/hg19_length.rda
- hg38: data/hg38_length.rda

## Workflow Example

```r
library(TiNDA)
data(hg19_length)

# Generate test data
test_df <- generate_test_data(hg19_length, num_variants = 500)

# Run TiNDA
tinda_obj <- TiNDA(test_df, sample_name = "sample_1", data_source = "WGS")

# Visualize
canopy_clst_plot(tinda_obj)
tinda_clst_plot(tinda_obj)
tinda_linear_plot(tinda_obj)
tinda_summary_plot(tinda_obj)

# Export results
write_data(tinda_obj, "results.tsv")
```

## Testing

No unit tests currently exist in the repository. Test files should be placed in `tests/testthat/`.

## Notes

- Input data should be pre-filtered rare germline variants (not common SNPs or artifacts)
- For large datasets, consider using only exonic variants
- Column names in input must match exactly: CHR, POS, Control_ALT_DP, Control_DP, Tumor_ALT_DP, Tumor_DP
