# ==============================================================================
# Lib_Plots.R
# requires following libraries
# - Lib_Analysis_Inversion.R
# ==============================================================================
# PROGRAMMERS:
# jb.feret@teledetection.fr
# jean-malo.perrot@inrae.fr
# ==============================================================================
# This Library includes functions aiming at producing figures to analyze
# performances of PROSPECT inversion
# ==============================================================================

#' This function produces a scatterplot comparing estimated value with measured
#' value, for one or multiple categories
#'
#' @param target numeric. Target values
#' @param estimate numeric. Estimated value to be compared to target values
#' @param Colors character. Vector of colornames corresponding to categories
#' @param Labs character. to be displayed as XY axes labels
#' @param fileName character. full path and name of the image to be saved
#' @param categories character. how to group inoput data. set to FALSE if no group
#' @param MinMaxAxis numeric. min and max values of XY axes
#'
#' @return none
#' @import ggplot2
#' @importFrom tools file_ext
#' @export

scatter_inversion <- function(target, estimate, Colors, Labs, fileName,
                              categories=FALSE,
                              MinMaxAxis = FALSE,
                              PlotStats = FALSE){
  
  # produce dataframe from input data
  input_df <- data.frame('target' = target,
                         'estimate' = estimate,
                         'categories' = categories)
  minmax <- c(min(c(target,estimate),na.rm = T),
              max(c(target,estimate),na.rm = T))
  
  if (isFALSE(MinMaxAxis[1])){
    MinMaxAxis <- minmax
  }
  # produce main plot
  Splot <- ggplot(input_df, aes(x=target, y=estimate, group=categories)) +
    geom_point(aes(pch = categories, color = categories,size=0.15, stroke = 0.25), size = 1.5) +
    xlim(MinMaxAxis[1], MinMaxAxis[2]) +
    ylim(MinMaxAxis[1], MinMaxAxis[2]) +
    scale_color_manual(values=Colors) +
    scale_shape_manual(values=c(19,19)) +
    labs(x=Labs[1],y=Labs[2]) +
    theme_bw() +
    theme(aspect.ratio=1,
          legend.position= "none")#,"bottom",
  #legend.title = element_text(color = "white"),
  #       legend.text = element_text(size = 5),
  #       axis.text = element_text(size=10),
  #       axis.title.x = element_text(size=10, face="bold"),
  #       axis.title.y = element_text(size=10, face="bold")) +
  # guides(fill=guide_legend(nrow = 2),size = "none")
  
  # Add 1:1 line
  Splot <- Splot + geom_abline(slope = 1, intercept = 0,linetype='dashed',size=0.75)
  
  # include statistics in the figure
  if (PlotStats == TRUE){
    nbCases <- length(unique(input_df$categories))
    case <- 0
    # compute statistics to compare performances
    Stats <- get_performances_inversion(target = input_df$target,
                                        estimate = input_df$estimate,
                                        categories = input_df$categories)
    for (conf in unique(1)){#BP_df[[BPvar]]$SpectralConf)){
      case <- case+1
      my.label.expr = paste0("italic(R)^ 2 == ", as.numeric(format(round(Stats$R2[case], 2), nsmall = 2)),sep='')
      Splot <- Splot + annotate(geom = "text",
                                x = MinMaxAxis[1] +0.02*(MinMaxAxis[2]-MinMaxAxis[1]),
                                y = MinMaxAxis[1] +(1.05-(0.07*case))*(MinMaxAxis[2]-MinMaxAxis[1]),
                                label = my.label.expr, parse = TRUE, hjust = 0,size=3, col="black")#Colors[case], fontface="bold")
      
      my.label.expr = paste0("italic(NRMSE) == ", as.numeric(format(round(Stats$NRMSE[case], 2), nsmall = 2)))
      Splot <- Splot + annotate(geom = "text",
                                x = MinMaxAxis[1] +0.02*(MinMaxAxis[2]-MinMaxAxis[1]),
                                y = MinMaxAxis[1] +(1.05-(0.17*case))*(MinMaxAxis[2]-MinMaxAxis[1]),
                                label = my.label.expr, parse = TRUE, hjust = 0, size=3, col= "black")#Colors[case], fontface="bold")
    }
  }
  # get iamge format
  ext_figure <- tools::file_ext(fileName)
  ggsave(filename = fileName, plot = Splot, device = ext_figure, path = NULL,
         scale = 1, width = 8, height = 8, units = "in", dpi = 600)
  return(invisible(Splot))
}
