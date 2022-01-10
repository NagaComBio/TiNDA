# Simulate test data for TiNDA

globalVariables('hg19_length')

#' Helper function to generate variant's depths
#' @param number_variants Number of variants to simulate
#' @param avg_depth Average depth/coverage for the variants
#' @param sd Standard deviation of the coverage
#' 
generate_depth = function(number_variants, avg_depth, sd = 5){
  return(abs(floor(rnorm(number_variants, mean = avg_depth, sd = sd))))
}

#' Simulate germline and somatic variants for a chromosome
#' @param chr chromosome name
#' @param chr_len chromosome length
#' @param avg_control_cov Average coverage of the control sample
#' @param avg_tumor_cov Average coverage of the tumor sample
#' @param avg_germline_MAF Average minor allele frequency (MAF) of 
#' the germline variants
#' @param avg_somatic_tumor_MAF Average MAF of the somatic variants 
#' in the tumor sample
#' @param avg_somatic_control_MAF Average MAF of the somatic variants in 
#' the control sample. This is the tumor in control (TiN) percentage, 
#' increase or decrease to simulate the TiN ratio.
#' @param num_variants Number of variants for a chromosome
#' @param per_somatic_variants Percentage of the variants to be somatic among 
#' the total number of 'num_variants'
#' @param ... ellipsis
#' 
simulate_variants <- function(chr, 
                              chr_len, 
                              avg_control_cov = 70, 
                              avg_tumor_cov = 70,
                              avg_germline_MAF = 0.50,
                              avg_somatic_tumor_MAF = 0.40,
                              avg_somatic_control_MAF = 0.03,
                              num_variants = 1000,
                              per_somatic_variants = 0.10, ...){
  
  total_variants = num_variants
  somatic_variants = num_variants * per_somatic_variants
  germline_variants = num_variants - somatic_variants
  
  # 
  tumor_germ_dp = generate_depth(germline_variants, avg_depth = avg_tumor_cov)
  tumor_germ_alt_dp = generate_depth(germline_variants, 
                                     avg_depth = avg_tumor_cov * avg_germline_MAF)
  
  tumor_som_dp = generate_depth(somatic_variants, avg_depth = avg_tumor_cov)
  tumor_som_alt_dp = generate_depth(somatic_variants, 
                                    avg_depth = avg_tumor_cov * avg_somatic_tumor_MAF)
  
  control_germ_dp = generate_depth(germline_variants, avg_depth = avg_control_cov)
  control_germ_alt_dp = generate_depth(germline_variants, 
                                       avg_depth = avg_control_cov * avg_germline_MAF)
  
  control_som_dp = generate_depth(somatic_variants, avg_depth = avg_control_cov)
  control_som_alt_dp = generate_depth(somatic_variants, 
                                      avg_depth = avg_control_cov * avg_somatic_control_MAF, 
                                      sd = 2)
  # For a chromosome
  test_data_chr <- tibble(CHR = rep(chr, total_variants),
                          POS = sample.int(chr_len, total_variants),
                          Control_ALT_DP = c(control_germ_alt_dp, control_som_alt_dp),
                          Control_DP = c(control_germ_dp, control_som_dp),
                          Tumor_ALT_DP = c(tumor_germ_alt_dp, tumor_som_alt_dp),
                          Tumor_DP = c(tumor_germ_dp, tumor_som_dp)
                          )
  
  return(test_data_chr)
}

#' Data format for TiNDA input
#' 
#' The input format is VCF like data frame object containing information 
#' from both tumor and normal samples. It is recommended to filter out 
#' common SNPs using gnomAD or other population databases and 
#' only use the rare germline variants. 
#' 
#' For user-defined input data, please use the same column names
#' as described below.
#' 
#' @format A data frame of rare germline variants
#' \itemize{
#'   \item CHR - chromosome name
#'   \item POS - variant position
#'   \item Control_ALT_DP - Read depth of the variant's alternate allele in the control sample
#'   \item Control_DP - Total read depth of the variant in the control sample
#'   \item Tumor_ALT_DP - Read depth of the variant's alternate allele in the tumor sample
#'   \item Tumor_DP - Total read_depth of the variant in the tumor sample
#' }
#' \tabular{rrrrrr}{
#' CHR \tab POS \tab Control_ALT_DP \tab Control_DP \tab Tumor_ALT_DP \tab Tumor_DP\cr
#' 1 \tab 1039001 \tab 20 \tab 40 \tab 23 \tab 46\cr
#' 1 \tab 2123023 \tab 12 \tab 32 \tab 14 \tab 23\cr
#' 1 \tab 3343543 \tab 23 \tab 56 \tab 34 \tab 67
#' }
#' 
#' @name Simulate_test_data
#' @docType data
#' 
#' @param chr_length_table A data frame containing chromosome names and its length
#' @param ... ellipsis
#' 
#' Simulate the test data generation for all the chromosomes
#' @export
generate_test_data <- function(chr_length_table, ...){
  chrs <- c(1:22, 'X', 'Y')
  test_data_all_chr <- lapply(chrs, function(x, ...){
    x_length = chr_length_table[chr_length_table$CHR == x,]$Length
    simulate_variants(x, x_length, ...)
  })
  
  test_data_all_chr <- bind_rows(test_data_all_chr)
  
  return(test_data_all_chr)
}

#' Chromosome length
#' 
#' Length of chromosomes for the linear TiNDA plots, default hg19. For used defined genomes, use the format below.
#' 
#' @usage data(hg19_length)
#' 
#' @name hg19_length
#' 
#' @format A data frame with chromosome name and its length
#' \itemize{
#'   \item CHR - chromosome name
#'   \item Length - Length of the chromosome
#' }
NULL
