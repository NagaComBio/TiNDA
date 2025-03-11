#' TiNDA function 
#' 
#' The main TiNDA function to rescue somatic variants from contaminated samples. 
#' Uses the EM-algorithm implemented in the Canopy R package.
#' It is recommend to use user's own local control database, 
#' MAF generated from local samples analyzed by the same pipeline, 
#' to filter our technical artifacts specific to their pipeline.
#' 
#' 
#' @param tbl A data frame object with rare germline variants with the raw coverage information.
#' @param sample_name Name of the sample for the title, Default: 'pid_1'
#' @param data_source WGS or WES, Default: 'WGS'
#' @param max_control_af Maximum control variant allele frequency (VAF), Default: 0.45
#' @param min_tumor_af Minimum tumor variant allele frequency (VAF), Default: 0.01
#' @param min_clst_members Minimum number of members in the cluster to be below the max_control_af and above the min_tumor_af, Default: 0.85
#' @param num_run For canopy cluster function, "number of EM runs for estimation for each specific number of clusters (to avoid EM being stuck in local optima)", Default: 1
#' @param min_control_af_chip Minimum control variant allele frequency (VAF) for CHIP cluster, Default: 0.02
#' @param max_control_af_chip Maximum control variant allele frequency (VAF) for CHIP cluster, Default: 0.35
#' @param max_tumor_af_chip Maximum tumor variant allele frequency (VAF) for CHIP cluster, Default: 0.25
#' @param find_chip Find CHIP clusters, Default: TRUE
#' @param verbose Print the potential clusters, Default: FALSE
#' @param ... ellipsis
#' 
#' @examples
#' data(hg19_length)
#' vcf_like_df = TiNDA::generate_test_data(hg19_length)
#' tinda_test_object <- TiNDA(vcf_like_df, sample_name = "sample_3", data_type = "WGS")
#' 
#' @importFrom stats median quantile rnorm
#' @importFrom methods is
#' @importFrom utils data
#' @importFrom magrittr "%>%"
#' @importFrom Canopy canopy.cluster
#' @import dplyr
#' @import purrr
#' 
#' @export
TiNDA <- function(tbl, 
                  sample_name = "pid_1",
                  data_source = "WGS",
                  max_control_af = 0.25,
                  min_tumor_af = 0.01,
                  min_clst_members = 0.85,
                  min_control_af_chip = 0.02,
                  max_control_af_chip = 0.40,
                  max_tumor_af_chip = 0.25,
                  num_run = 1,
                  find_chip = TRUE,
                  verbose = FALSE, ...) {
  
  if(data_source == "WGS") {
    mu.init <- cbind(c(0.5, 0.95, 0.50, 0.50, 0.02, 0.02, 0.02, 0.25, 0.15, 0.35), 
                     c(0.5, 0.95, 0.15, 0.85, 0.20, 0.50, 0.85, 0.02, 0.02, 0.02))
    numberCluster <- 10
  } else if(data_source == "WES") {
    mu.init <- cbind(c(0.5, 0.95, 0.50, 0.50, 0.02, 0.02, 0.02, 0.25, 0.15, 0.35), 
                     c(0.5, 0.95, 0.15, 0.85, 0.20, 0.50, 0.85, 0.02, 0.02, 0.02))
    numberCluster <- 10
  }
  
  # Test tbl 
  expected_col_names <- c('CHR', 'POS', 'Control_ALT_DP', 'Control_DP',
                          'Tumor_ALT_DP', 'Tumor_DP')
  if (!purrr::is_empty(setdiff(expected_col_names, colnames(tbl)))){
    cat("Expected column names didn't appear in the input table\n")
    cat("Expected:", expected_col_names, "\n")
    stop()
  }
  
  
  cat("Found ", dim(tbl)[1], " variants from the input data\n")
  new_tbl <- tbl %>% 
    filter(.data$Control_ALT_DP < .data$Control_DP,
           .data$Tumor_ALT_DP < .data$Tumor_DP) -> new_tbl
  
  if(dim(tbl)[1] == dim(new_tbl)[1]){
    rm(new_tbl)
  } else {
    num_removed <- dim(tbl)[1] - dim(new_tbl)[1]
    cat("Removed ", num_removed)
    cat(" variants with ALT read depth more than total read depth\n")
    tbl = new_tbl
  }
  # Variant AF
  tbl %>% 
    mutate(Control_AF = .data$Control_ALT_DP/ .data$Control_DP,
           Tumor_AF = .data$Tumor_ALT_DP/ .data$Tumor_DP) -> tbl
  
  # Running Canopy -------------------------------------------------------------
  R <-as.matrix(tbl[,c('Control_ALT_DP', 'Tumor_ALT_DP')])
  X <-as.matrix(tbl[,c('Control_DP', 'Tumor_DP')])

  # Canopy run and assigning centers
  try(
      canopy.clust <- canopy.cluster(R, X, 
                                     num_cluster = numberCluster, 
                                     num_run = num_run, Mu.init = mu.init)
    )
  if(is.null(canopy.clust)){
    stop("Failed canopy run\n", call. = FALSE)
  }
  tbl$canopyCluster<-canopy.clust$sna_cluster

  # Select the potential TiN clusters ------------------------------------------
  potential_somatic_clst_per <- tbl %>%  
    mutate(squareRescue = .data$Control_AF < max_control_af & 
                          .data$Tumor_AF > min_tumor_af) %>% 
    mutate(diagonalRescue = .data$Control_AF < .data$Tumor_AF) %>% 
    group_by(.data$canopyCluster) %>%
    dplyr::summarise(prop1 = mean(.data$squareRescue == TRUE), 
                     prop2 = mean(.data$diagonalRescue==TRUE))

  potential_somatic_clst <- potential_somatic_clst_per %>% 
    filter(.data$prop1 >  min_clst_members & 
             .data$prop2 > min_clst_members) %>% 
    dplyr::select(.data$canopyCluster) %>% 
    collect() %>% pull(.data$canopyCluster)
  
    # in verbose mode, print the potential clusters
  if(verbose){
    cat("Potential somatic clusters\n")
    print(potential_somatic_clst_per)
    cat("Selected somatic clusters\n")
    print(potential_somatic_clst)
  }
  
  # Clonal hematopoiesis clusters ------------------------------------------
  if(!find_chip){
    potential_chip_clst <- c()
  } else {
    potential_chip_clst_per <- tbl %>%
      mutate(squareRescue = .data$Control_AF > min_control_af_chip & 
                            .data$Control_AF < max_control_af_chip & 
                            .data$Tumor_AF < max_tumor_af_chip) %>% 
      mutate(diagonalRescue = .data$Control_AF > .data$Tumor_AF) %>%
      group_by(canopyCluster) %>%
      dplyr::summarise(prop1 = mean(squareRescue == T),
                      prop2 = mean(diagonalRescue==T))

    potential_chip_clst <- potential_chip_clst_per %>%
      filter(.data$prop1 > min_clst_members & 
              .data$prop2 > min_clst_members) %>%
      dplyr::select(.data$canopyCluster) %>%
      collect() %>% pull(.data$canopyCluster)
    
    # in verbose mode, print the potential clusters
    if(verbose){
      cat("Potential CHIP clusters\n")
      print(potential_chip_clst_per)
      cat("Selected CHIP clusters\n")
      print(potential_chip_clst)
    }
  }
  
  # Potential germline clusters -----------------------------------------------
  potential_germline_clst <- tbl %>%
    filter(!.data$canopyCluster %in% c(potential_somatic_clst, potential_chip_clst)) %>%
    dplyr::select(.data$canopyCluster) %>%
    collect() %>% pull(.data$canopyCluster)
  
  # Removing unusual clusters ---------------------------------------------------
  # Usually occurs due to presence of variants in-between somatic and true germline clusters
  potential_somatic_clst_copy <- potential_somatic_clst

  for (cl in potential_somatic_clst_copy) {
    # Get the cluster size
    cluster_size <- nrow(tbl[tbl$canopyCluster == cl, ])
    
    # Find the control line with the maximum Control_AF in the cluster
    control_max_line <- tbl[tbl$canopyCluster == cl, ][which.max(tbl[tbl$canopyCluster == cl, ]$Control_AF), ]
    # Count the number of mutations in other clusters that meet the criteria
    count <- nrow(tbl[tbl$Control_AF <= control_max_line$Control_AF & 
                        tbl$Tumor_AF > control_max_line$Tumor_AF & 
                      !(tbl$canopyCluster %in% potential_somatic_clst), ])
  
    
    # Remove the cluster from the list if the count is greater than a quarter of the cluster size
    if (count > cluster_size / 2) {
      potential_somatic_clst <- potential_somatic_clst[!(potential_somatic_clst %in% cl)]
    }
  }

  # Rescuing the TiN homozygous cluster at the top-left corner -----------------
  if(length(potential_somatic_clst) != 0) {
    # Finding the top most cluster in TiN rescue region
    tbl %>% group_by(.data$canopyCluster) %>% 
        summarise(max_C = max(.data$Control_AF), 
                  med_T = quantile(.data$Tumor_AF, 0.80)) %>% 
        filter(.data$canopyCluster %in% potential_somatic_clst ) %>%
        filter(.data$med_T == max(.data$med_T)) -> homo.threshold
  
    # Counting the total number of TiN variants if we rescue the homozygous variants we well
    tbl %>% filter((.data$Tumor_AF > homo.threshold$med_T & 
                      .data$Control_AF < homo.threshold$max_C) | 
                     .data$canopyCluster %in% potential_somatic_clst) %>% 
      nrow() -> TiN.homo.updated
    
    # Only the variants from Tin Clusters
    tbl %>% filter(.data$canopyCluster %in% potential_somatic_clst) %>% 
      nrow() -> TiN.not.homo.updated
  
    # The Rescued TiN should only add 10% more variants
    if(TiN.homo.updated <= (TiN.not.homo.updated + (TiN.not.homo.updated/10))) {
      tbl %>% 
        mutate(TiN_Class = ifelse((.data$Tumor_AF > homo.threshold$med_T &
                                     .data$Control_AF < homo.threshold$max_C) | 
                                    .data$Control_AF == 0 | 
                                  .data$canopyCluster %in% potential_somatic_clst, 
                                  "Somatic_Rescue", "Germline")) -> tbl
    } else {
      tbl %>% mutate(TiN_Class = ifelse(.data$Control_AF == 0 | 
                                          .data$canopyCluster %in% potential_somatic_clst, 
                                        "Somatic_Rescue", "Germline")) -> tbl
    } 
  } else {
    tbl %>% mutate(TiN_Class = ifelse(.data$Control_AF == 0, 
                                      "Somatic_Rescue", "Germline")) -> tbl
  }

  # Annotating CHiP Clusters
  if(length(potential_chip_clst) != 0){
    tbl %>%
    mutate(
      TiN_Class = ifelse(.data$canopyCluster %in% potential_chip_clst & !grepl("Somatic_", TiN_Class), "CHIP", TiN_Class)
    ) -> tbl
}

  tinda_object <- list(data=tbl, 
                       min_tumor_af = min_tumor_af,
                       max_control_af = max_control_af,
                       number_cluster = numberCluster,
                       sample_name = sample_name)
  
  class(tinda_object) <-'TiNDA'
  
  return(tinda_object)
}
