  <!-- badges: start -->
  [![R build status](https://github.com/NagaComBio/TiNDA/workflows/R-CMD-check/badge.svg)](https://github.com/NagaComBio/TiNDA/actions)
  [![CRAN version](https://www.r-pkg.org/badges/version/TiNDA)](https://cran.r-project.org/package=TiNDA)
  <!-- badges: end -->

# TiNDA
## Tumor in Normal Detection Analysis

## Overview

This is an R package to rescue somatic variants called as germline due to tumor DNA contamination in the patient's blood/control sample.

TiNDA makes use of the [Canopy's](https://github.com/yuchaojiang/Canopy) EM-cluster function to partition the variants into different clusters. And uses the following assumptions to define these clusters into somatic and germline.

Based on the following assumptions:

1. The variant allele frequency (VAF) of somatic variants in tumor samples will be higher than contaminated somatic variants in the control sample.
2. The contamination exceeding a certain threshold (`max_control_af: 0.25`) will be difficult to separate from the germline VAF.

An area of interest (AOI) is defined in the control vs tumor VAF 2D space. Clusters with a majority (`min_clst_members: 0.85`) of its members within this AOI are defined as 'omatic rescue'.

## Area of Interest
In the tumor VAF vs control VAF, the AOI for somatic and ChiP variants are defined in the following image. The "golden" polygon defines the somatic region, and the "red" polygon defines the ChiP region, with the rest of the areas defining germline variants.

![AOI](man/figures/polygon_aoi.png)


## Key Features
- **Rescuing Misclassified Variants**: TiNDA rescues somatic variants that are misclassified due to tumor-in-normal contamination.
- **Detecting CHiP Clusters**: TiNDA identifies CHiP clusters by distinguishing germline variants from genuine somatic mutations in blood.
- **Visualization**: TiNDA provides visualization tools to help users assess quality of the clustering.


## Installation

Install directly from the GitHub

```
devtools::install_github("nagacombio/tinda")
```

## Usage

#### Workflow
The TiNDA input consists of read counts for rare and private variants, including both germline and somatic variants. These variants should be identified through the joint analysis of tumor and control samples, and they must be filtered to remove common SNPs and technical artifacts. If the dataset is still too large and to expedite clustering and plotting, consider using only exonic variants.
 
An ideal workflow with TiNDA:

![TiNDA workflow](man/figures/tinda_flow.png)

#### Input data format

The input data for TiNDA is a data frame containing the following information/columns, 

  * **CHR** - Chromosome name
  * **POS** - Variant position
  * **Control_ALT_DP** - Read depth of the variant's alternate allele in the control sample
  * **Control_DP** - Total read depth of the variant in the control sample
  * **Tumor_ALT_DP** - Read depth of the variant's alternate allele in the tumor sample
  * **Tumor_DP** - Total read_depth of the variant in the tumor sample

**Note:** Keep the column names in the input table.

An example table,

|CHR| POS | Control_ALT_DP | Control_DP | Tumor_ALT_DP | Tumor_DP
|--|--:|--:|--:|--:|--:
 1 | 1039001 | 20 | 40 | 23 | 46
 1 | 2123023 | 12 | 32 | 14 | 23
 1 | 3343543 | 23 | 56 | 34 | 67

#### Example TiNDA analysis
```{r}
# Generate data to test the package
library(TiNDA)
data(hg19_length)
test_df <- generate_test_data(hg19_length, num_variants = 500)
```

#### Run the TiNDA function

```{r}
# Check the documentation for the paramaters
tinda_object <- TiNDA(test_df)
```

#### Plotting the results
```{r}
# Plot the results of the canopy cluster analysis
canopy_clst_plot(tinda_object)
```
![canopy_clst_plot](man/figures/canopy_clst_plot.png)

```{r}
# Plot the TiNDA cluster assignment
tinda_clst_plot(tinda_object)
```
![tinda_clst_plot](man/figures/tinda_clst_plot.png)

```{r}
# Plot the linear plot of the TiNDA results
tinda_linear_plot(tinda_object)
```
![tinda_linear_plot](man/figures/tinda_linear_plot.png)

```{r}
# Plot the summary of the TiNDA results - includes canopy clusters, TiNDA cluster assignment and linear plots
tinda_summary_plot(tinda_object)
```
![tinda_summary_plot](man/figures/tinda_summary_plot.png)

---

## Changelog

### [1.2.0] - 2026-04-29

#### Added
- `get_tinda_params()`: Function to retrieve default parameters for WGS/WES analysis
- `run_pipeline()`: Convenience function for file-based workflow (input → analysis → output)
- `print.TiNDA()`: S3 print method for TiNDA objects
- `summary.TiNDA()`: S3 summary method providing detailed classification statistics
- Test infrastructure with testthat and 10+ basic tests
- `inst/CITATION`: Citation file for publication credit

#### Changed
- Enhanced TiNDA object structure with `classification_summary` and `parameters` fields
- Improved error messages for input validation (missing columns, invalid data types)

#### Fixed
- Documentation: `max_control_af` default value (0.45 → 0.25)
- Documentation: `max_control_af_chip` default value (0.35 → 0.40)
- Example code: `data_type` → `data_source` parameter name in all examples
- Input validation: silent failures now throw informative errors

#### Security
- Added comprehensive input validation to prevent invalid data processing

---

### [1.1.0] - Previous Release

#### Added
- Initial release with core TiNDA analysis functionality
- Canopy-based EM clustering for variant classification
- Visualization functions: `canopy_clst_plot()`, `tinda_clst_plot()`, `tinda_linear_plot()`, `tinda_summary_plot()`
- Data simulation functions: `generate_test_data()`, `simulate_variants()`, `generate_depth()`
- Output function: `write_data()`
- Support for hg19 and hg38 reference genomes

---

## Versioning

This project follows [Semantic Versioning](https://semver.org/):

- **MAJOR** version: Incompatible API changes
- **MINOR** version: Backward-compatible new functionality
- **PATCH** version: Backward-compatible bug fixes

For the latest version information, check the GitHub releases page.