assert_class <- function(obj){
  if(!is(obj, "TiNDA")) stop("The input object is not of class 'TiNDA'", call. = F)
}


#' Plotting the Canopy identified clusters
#' 
#' Using ggplot2, plots VAF from tumor and control, and colors the based on canopy cluster
#' 
#' @param tinda_object Object returned by TiNDA function
#' 
#' @examples  
#' tinda <- TiNDA(coverage_object, sample_name = "sample_3", data_type = "WGS")
#' canopy_clst_plot(tinda)
#'
#' @import ggplot2
#' 
#' @export

canopy_clst_plot <- function(tinda_object, ...){
  
  assert_class(tinda_object)
  
  tbl <- tinda_object$data
  max_control_af <- tinda_object$max_control_af
  min_tumor_af <- tinda_object$min_tumor_af
  
  poly.df <- data.frame(x=c(0, 0, max_control_af, max_control_af), 
                        y=c(min_tumor_af, 1, 1, max_control_af))
  ggplot() + 
    geom_point(data=tbl, aes(Control_AF, Tumor_AF, color=factor(canopyCluster)), alpha=0.5) + 
    geom_polygon(data=poly.df, aes(x, y), alpha=0.2, fill="gold") + 
    theme_bw() + 
    theme(text = element_text(size=15), legend.position="bottom") + 
    xlab("Control VAF") + 
    ylab("Tumor VAF") + 
    xlim(0,1) + ylim(0,1) +
    guides(color=guide_legend("Canopy clusters")) + 
    ggtitle(paste0("Clusters from Canopy"))
}

#' Plotting the  TiNDA indentified TiN clusters
#' 
#' Plotting the somatic resuced and germline clusters as classified by TiNDA
#' 
#' @param tinda_object Object returned by TiNDA function
#' 
#' @examples 
#' tinda <- TiNDA(coverage_object, sample_name = "sample_3", data_type = "WGS")
#' tinda_clst_plot(tinda)
#' 
#' @import ggplot2
#' 
#' @export


tinda_clst_plot <- function(tinda_object, ...){
  
  assert_class(tinda_object)
  
  tbl <- tinda_object$data
  max_control_af <- tinda_object$max_control_af
  min_tumor_af <- tinda_object$min_tumor_af
  
  poly.df <- data.frame(x=c(0, 0, max_control_af, max_control_af), 
                        y=c(min_tumor_af, 1, 1, max_control_af))
  
  ggplot() + 
    geom_point(data=tbl, aes(Control_AF, Tumor_AF, color=TiN_Class), alpha=0.3) + 
    geom_polygon(data=poly.df, aes(x, y), alpha=0.2, fill="gold") + 
    theme_bw() + 
    theme(text = element_text(size=15), legend.position="bottom") + 
    xlab("Control VAF") + 
    ylab("Tumor VAF") + 
    xlim(0,1) + ylim(0,1) +
    guides(color=guide_legend("TiN clusters")) + 
    ggtitle(paste0("TiN clusters"))
}


#' Plotting linear genome
#' 
#' Function to plot choromosomes linerlly from TiNDA object
#' 
#' @param tinda_object  Object returned by TiNDA function
#' @param sample_af sample AF column in the TiNDA data, Default: 'Tumor_AF'
#' @param chr.length Data frame object of length of the chrs, if null hg19 chr length is used, Default: NULL'  
#' @param colorCol Name of column with the final TiNDA classification, categorical variable, Default: 'TiN_Class'
#' 
#' @examples 
#' tinda <- TiNDA(coverage_object, sample_name = "sample_3", data_type = "WGS")
#' plotGenome_ggplot(tinda)
#' plotGenome_ggplot(tinda, sample_af = 'Control_AF')
#' 
#' @import ggplot2
#' 
#' @export

tinda_linear_plot <- function(tinda_object, 
                              sample_af = 'Tumor_AF', 
                              chr.length = NULL, 
                              colorCol = 'TiN_Class', ...) {
  assert_class(tinda_object)
  
  if(is.null(chr.length)){
    data("hg19_length")
    chr.length <- hg19_length
    chr.length$newShift<-c(0, chr.length$Length[-length(chr.length$Length)])
    chr.length$cumShiftLength <- cumsum(as.numeric(chr.length$new))
    chr.length$cumLength <- cumsum(as.numeric(chr.length$Length))
    chrTotalLength <- sum(as.numeric(chr.length$Length))
  }
  
  tbl <- tinda_object$data  
  tbl %>% 
    mutate(CHR=as.factor(CHR)) %>% 
    left_join(chr.length, by="CHR") %>% 
    mutate(cumPOS = POS + cumShiftLength) %>% 
    ggplot() + 
    geom_point(aes_string(x="cumPOS", y=sample_af, color=colorCol), shape=124) + 
    theme_bw() + 
    scale_x_continuous(breaks = chr.length$cumLength - (chr.length$Length/2),
                       labels = chr.length$CHR,
                       minor_breaks = chr.length$cumShiftLength,
                       expand = c(0, 0)) + 
    theme(axis.title.x = element_blank(), text = element_text(size=15),
          panel.grid.major.x = element_line(color="lightgrey", linetype = 0),
          panel.grid.minor.x = element_line(color="grey")) 
}

#' TiNDA summary plot
#' 
#' A summary plot combines differnt TiNDA plots and TiNDA results
#' 
#' @param tinda_object Object returned by TiNDA function
#' 
#' @examples 
#' 
#' @import grid
#' @import gridExtra
#' 
#' @export


tinda_summary_plot <- function(tinda_object, ...){
  
  assert_class(tinda_object)
  
  tbl <- tinda_object$data
  sample_name <- tinda_object$sample_name
    
  tinda_summary_info <- tbl %>% 
    group_by(TiN_Class) %>% 
    summarise(Count =n(), Median_Control_VAF = formatC(median(Control_AF), digits=5, format="f"),
              Median_Tumor_VAF = formatC(median(Tumor_AF), digits=5, format="f"))
  
  TableTheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 1, hjust=1, x=0.95)),
    colhead = list(fg_params=list(cex = 1, hjust=1, x=0.95)),
    rowhead = list(fg_params=list(cex = 1, hjust=1, x=0.95)))
  
  TableAnn <-tableGrob(tinda_summary_info, rows = c(), theme=TableTheme)
  
  PlotLayout <-rbind(c(1,2,3),
                     c(1,2,3),
                     c(4,4,4),
                     c(5,5,5))
  
  p1 <- canopy_clst_plot(tinda_object)
  p2 <- tinda_clst_plot(tinda_object)
  p3 <- tinda_linear_plot(tinda_object)
  p4 <- tinda_linear_plot(tinda_object, sample_af = 'Control_AF')
  
  grid.arrange(p1, p2, TableAnn, p3, p4, 
               layout_matrix = PlotLayout,
               top=textGrob(paste0('Tumor in Normal Detection Analysis - TiNDA : ', sample_name), 
                            gp=gpar(cex=2)))
}
