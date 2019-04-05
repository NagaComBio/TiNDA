#' TiNDA function 
#' 
#' The main TiNDA function to rescue somatic variants from contaminated samples. 
#' Uses the EM-algorithm implemented in the Conopy R package.
#' It is recommend to use user's own local control database, 
#' MAF generated from local samples analyzed by the same pipeline, 
#' to filter our techincal artifacts specific to their pipeline.
#' 
#' 
#' @param tbl A data frame object with rare germline variants with raw coverage infomation.
#' @param sample_name Name of the sample for the title, Default: 'pid_1'
#' @param data_source WGS or WES, Default: 'WGS'
#' @param max_control_af Maximum control variant allele frequency (VAF), Default: 0.45
#' @param min_tumor_af Minimum tumor varint allele frequency (VAF), Default: 0.01
#' @param min_clst_members Minimum number of members in the cluster to be below the max_control_af and above the min_tumor_af, Default: 0.85
#' @param num_run For canopy cluster function, "number of EM runs for estimation for each specific number of clusters (to avoid EM being stuck in local optima)", Default: 1
#' 
#' @examples 
#' TiNDA(vcf_like_df, sample_name = "sample_1", data_type = "WES", max_control_af = 0.35)
#' TiNDA(vcf_like_df, sample_name = "sample_3", data_type = "WGS")
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Canopy  canopy.cluster
#' @import dplyr
#' 
#' @export
TiNDA <- function(tbl, 
                  sample_name = "pid_1",
                  data_source = "WGS",
                  max_control_af = 0.45,
                  min_tumor_af = 0.01,
                  min_clst_members = 0.85,
                  num_run = 1, ...) {
  
  if(data_source == "WGS") {
    mu.init <- cbind(c(0.5, 0.95, 0.5,  0.5,  0.5,  0.5,  0.02, 0.02, 0.02, 0.02, 0.10), 
                    c(0.5, 0.95, 0.25, 0.75, 0.95, 0.05, 0.30, 0.5,  0.95, 0.10, 0.10))
    numberCluster <- 11
  } else if(data_source == "WES") {
    mu.init <- cbind(c(0.5, 0.02, 0.02, 0.1, 0.02), 
                    c(0.5, 0.30, 0.5,  0.1, 0.10))
    numberCluster <- 5
  }
  
  # Test tbl 
  cat("Found ", dim(tbl)[1], " variants from the input data\n")
  new_tbl <- tbl[tbl[,3] <= tbl[,4],]
  new_tbl <- new_tbl[new_tbl[,5] <= new_tbl[,6],]
  if(dim(tbl)[1] == dim(new_tbl)[1]){
    rm(new_tbl)
  } else {
    num_removed <- dim(tbl)[1] - dim(new_tbl)[1]
    cat(num_removed, "variants with alternate read depth more than total read depth were removed\n")
    tbl = new_tbl
  }
  # Variant AF 
  tbl$Control_AF = tbl[,3]/tbl[,4]
  tbl$Tumor_AF   = tbl[,5]/tbl[,6]
  
  # Running Canopy ---------------------------------------------------------------
  R <-as.matrix(tbl[,c(5,3)]) # Tumor read depth
  X <-as.matrix(tbl[,c(6,4)]) # Control read depth

  print(head(R))
  print(head(X))
  # Canopy run and assigning centers
  try(
      canopy.clust <- canopy.cluster(R, X, num_cluster = numberCluster, 
                               num_run = num_run, Mu.init = mu.init)
    )
  if(is.null(canopy.clust)){
    stop("Failed canopy run\n", call. = F)
  }
  tbl$canopyCluster<-canopy.clust$sna_cluster

  # Select the potential TiN clusters -------------------------------------------
  potential_somatic_clst <- tbl %>%  
    mutate(squareRescue = Control_AF < max_control_af & 
             Tumor_AF > min_tumor_af) %>% 
    mutate(diagonalRescue = Control_AF < Tumor_AF) %>% 
    group_by(canopyCluster) %>%
    dplyr::summarise(prop1 = mean(squareRescue == T), 
                     prop2 = mean(diagonalRescue==T)) %>% 
    filter(prop1 >  min_clst_members & 
             prop2 > min_clst_members) %>% 
    dplyr::select(canopyCluster) %>% 
    collect() %>% .[[1]]

  # Removing unusual clusters ---------------------------------------------------
  # Usually occurs due to presence of variants in-between somatic and true germline clusters
  for(cl in potential_somatic_clst){
    cluster.size <- nrow(tbl[tbl$canopyCluster == cl,])
    control.max.line <- tbl[tbl$canopyCluster == cl,][which.max(tbl[tbl$canopyCluster == cl,]$Control_AF),]
  
    count <- nrow(tbl[tbl$Control_AF <=  control.max.line$Control_AF & 
                        tbl$Tumor_AF > control.max.line$Tumor_AF & 
                        !(tbl$canopyCluster %in% potential_somatic_clst),])
  
    if(count > cluster.size/4){
      potential_somatic_clst <- potential_somatic_clst[!potential_somatic_clst == cl]
    }
  }

  # Rescuing the TiN homozygous cluster at the top-left corner ----------------------
  if(length(potential_somatic_clst) != 0) {
    # Finding the top most cluster in TiN rescue region
    tbl %>% group_by(canopyCluster) %>% 
        summarise(max_C = max(Control_AF), 
                  med_T = quantile(Tumor_AF, 0.80)) %>% 
        filter(canopyCluster %in% potential_somatic_clst ) %>%
        filter(med_T == max(med_T)) -> homo.threshold
  
    # Couting the total number of TiN variants if we rescue the homozygous variants we well
    tbl %>% filter((Tumor_AF > homo.threshold$med_T & Control_AF < homo.threshold$max_C) | 
                     canopyCluster %in% potential_somatic_clst) %>% 
      nrow() -> TiN.homo.updated
    
    # Only the variants from Tin Clusters
    tbl %>% filter(canopyCluster %in% potential_somatic_clst) %>% 
      nrow() -> TiN.not.homo.updated
  
    # The Rescued TiN should only add 10% more variants
    if(TiN.homo.updated <= (TiN.not.homo.updated + (TiN.not.homo.updated/10))) {
      tbl %>% 
        mutate(TiN_Class = ifelse((Tumor_AF > homo.threshold$med_T & Control_AF < homo.threshold$max_C) | 
                                    Control_AF == 0 | 
                                  canopyCluster %in% potential_somatic_clst, "Somatic_Rescue", "Germline")) -> tbl
    } else {
      tbl %>% mutate(TiN_Class = ifelse(Control_AF == 0 | 
                                          canopyCluster %in% potential_somatic_clst, "Somatic_Rescue", "Germline")) -> tbl
    } 
  } else {
    tbl %>% mutate(TiN_Class = ifelse(Control_AF == 0, "Somatic_Rescue", "Germline")) -> tbl
  }

  tinda_object <- list(data=tbl, 
                       min_tumor_af = min_tumor_af,
                       max_control_af = max_control_af,
                       number_cluster = numberCluster,
                       sample_name = sample_name)
  
  class(tinda_object) <-'TiNDA'
  
  return(tinda_object)
}
