# ==============================================================================

# This library makes it possible to simulate an optical sensor for the 

# This library makes it possible to simulate an optical sensor for the

# physiological monitoring of leaves from the PROSPECT model. Each function is
# meant to be as generalized as possible
# Lib_Gaussian_filter.R
# ==============================================================================
# PROGRAMMERS: Jean-Baptiste FERET, Jean-Malo PERROT
#
# ==============================================================================
# This Library includes nmrse_func
# ==============================================================================

#' This function allows application of a gaussian filter on reflectance and
#' transmittance data
#'
#' @param obs list(num) : observed data.
#' @param pred list(num) : predicted data.
#' @param type numeric : type of normalization .
#'
#' @return nrmse : Normalized RMSE
#' @importFrom 
#' @export
#' @importFrom
#' @export
#' 
nrmse_func <-  function(obs, pred, type = "sd") {
  #from https://www.marinedatascience.co/blog/2019/01/07/normalizing-the-rmse/

  squared_sums <- sum((obs - pred)^2)
  mse <- squared_sums/length(obs)                                               #Mean squared error
  rmse <- sqrt(mse)                                                             #Root mean squared error
  if (type == "sd") nrmse <- rmse/sd(obs)                                       #Normalized RMSE by standard deviation
  if (type == "mean") nrmse <- rmse/mean(obs)                                   #Normalized RMSE by mean
  if (type == "maxmin") nrmse <- rmse/ (max(obs) - min(obs))                    #Normalized RMSE by max-min
  if (type == "iq") nrmse <- rmse/ (quantile(obs, 0.75) - quantile(obs, 0.25))  #Normalized RMSE by inter-quartile range
  
  return(nrmse)
}

