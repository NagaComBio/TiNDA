#' Data format for TiNDA input
#' 
#' The input format is VCF like data frame object containing information from both tumor and normal samples.
#' It is recommened to filter out common SNPs using gnomAD or other population databases and use only the rare germline variants
#' 
#' @format A data frame of rare germline variants
#' \itemize{
#'   \item CHR - chromsome name
#'   \item POS - variant position
#'   \item CONTROL_ALT_DP - Read depth of the variant's alternate allele in the control sample
#'   \item CONTROL_DP - Total read depth of the variant in the control sample
#'   \item TUMOR_ALT_DP - Read depth of the variant's alternate allele in the tumor sample
#'   \item TUMOR_DP - Total read_depth of the variant in the tumor sample
#' }
#' \tabular{rrrrrr}{
#' CHR \tab POS \tab CONTROL_ALT_DP \tab CONTROL_DP \tab TUMOR_ALT_DP \tab TUMOR_DP\cr
#' 1 \tab 1039001 \tab 20 \tab 40 \tab 23 \tab 46\cr
#' 1 \tab 2123023 \tab 12 \tab 32 \tab 14 \tab 23\cr
#' 1 \tab 3343543 \tab 23 \tab 56 \tab 34 \tab 67
#' }
#' 
#' @name vcf_like_df
#' @docType data
NULL

#' Chromosome length
#' 
#' Length of chromosomes for the linear TiNDA plots, default hg19. For used defined genomes, use the format below.
#' 
#' @usage data(hg19_length)
#' 
#' @name hg19_length
#' 
#' @format A data frame with chr name and its length
#' \itemize{
#'   \item CHR - chromosome name
#'   \item Length - Length of the chromosome
#' }
NULL