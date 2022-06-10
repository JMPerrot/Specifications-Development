# ==============================================================================
# Lib_Analysis_Inversion.R
# ==============================================================================
# PROGRAMMERS:
# jb.feret@teledetection.fr
# jean-malo.perrot@inrae.fr
# ==============================================================================
# This Library includes functions aiming at analyzing results of PROSPECT inversion
# ==============================================================================

#' This function computes various statistical indices to assess the performances
#' of the estimations obtained from a model compared to expected values
#'
#' @param target numeric. Target values
#' @param estimate numeric. Estimated value to be compared to target values
#' @param categories character. how to group inoput data. set to FALSE if no group
#'
#' @return Stats_df data.frame. Includes all statistical indicators computed for each category
#' @importFrom stats IQR
#' @importFrom Metrics rmse
#' @importFrom stats lm
#' @export



get_performances_inversion <- function(target, estimate, categories=FALSE){

  # produce dataframe from input data
  input_df <- data.frame('target' = target,
                         'estimate' = estimate,
                         'categories' = categories)

  # compute statistics for the resulting data frame
  R2 <- RMSE <- NRMSE <- list()
  # compute IQR for NRMSE
  # https://www.marinedatascience.co/blog/2019/01/07/normalizing-the-rmse/
  iqr <- stats::IQR(input_df$target) #interquartil ratio

  for (conf in unique(input_df$categories)){
    subset <- input_df %>% filter(input_df$categories==conf)
    # compute coefficient of determination
    linmod = lm(target ~ estimate, data=subset)
    R2[[conf]] <- summary(linmod)$r.squared
    # compute NRMSE
    RMSE[[conf]] <- Metrics::rmse(actual = subset$target, predicted = subset$estimate)
    NRMSE[[conf]] <- 100*RMSE[[conf]]/iqr
  }
  Stats_df <- data.frame('R2'=unlist(R2), 'RMSE'= unlist(RMSE), 'NRMSE'=unlist(NRMSE))
return(Stats_df)
}
