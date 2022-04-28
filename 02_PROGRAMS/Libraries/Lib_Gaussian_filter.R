# ==============================================================================
# This library makes it possible to simulate an optical sensor for the 
# physiological monitoring of leaves from the PROSPECT model. Each function is
# meant to be as generalized as possible
# Lib_Gaussian_filter.R
# ==============================================================================
# PROGRAMMERS: Jean-Baptiste FERET, Jean-Malo PERROT
#
# ==============================================================================
# This Library includes gaussian_filter
# ==============================================================================

#' This function allows application of a gaussian filter on reflectance and
#' transmittance data
#'
#' @param fwhm numeric : full width at half maximum.
#' @param wl list(num) : wavelength of measured data.
#' @param SRF list : Sensor Response Function.
#' @param Refl Large matrix(num) : list of reflectance data.
#' @param Tran Large matrix(num) : list of transmittance data.
#'
#' @return LRT_df data.frame : dataframe with wavelength, reflectance and transmittance
#' @importFrom prosail applySensorCharacteristics
#' @export


gaussian_filter <- function(SRF, lambda, Refl=NULL, Tran=NULL){
  # apply SRF to leaf optical properties from ANGERS dataset
  if (!is.null(Refl)){
    Refl_filter <- prosail::applySensorCharacteristics(wvl = lambda,
                                                       InRefl = Refl,
                                                       SRF = SRF)
  } else {
    Refl_filter <- NULL
  }
  if (!is.null(Tran)){
    Tran_filter <- prosail::applySensorCharacteristics(wvl = lambda,
                                                     InRefl = Tran,
                                                     SRF = SRF)
  } else {
    Tran_filter <- NULL
  }
  LRT_df <- list('wl'=SRF$Spectral_Bands,'Refl'=Refl_filter,'Tran'=Tran_filter)
  return(LRT_df)
}


#' This function performs PROSPECT inversion using updated spectral sampling and
#' gaussian filter defined by FHWM
#'
#' @param Sampling_in numeric. spectral sampling of input data
#' @param Refl numeric. Reflectance with Sampling_in
#' @param Tran numeric. Transmittance with Sampling_in
#' @param fwhm numeric. FHWM of output optical properties
#' @param minWL_OUT numeric. minimum spectral band for updated sampling
#' @param maxWL_OUT numeric. maximum spectral band for updated sampling
#'
#' @return
#' @importFrom prosail applySensorCharacteristics
#' @importFrom progress progressbar
#' @export

Invert_PROSPECT_GaussianFilter <- function(Sampling_in, Refl, Tran, fwhm, minWL_OUT, maxWL_OUT){

  minWL <- min(Sampling_in)#min(lambda)
  maxWL <- max(Sampling_in)#max(lambda)
  # define spectral sampling based on minWL, minWL_OUT, FWHM
  wl <- seq(max(minWL_OUT,minWL+(2*fwhm)),
            min(maxWL_OUT,maxWL-(2*fwhm)),
            by = fwhm)                           #wavelength range
  SpectralProps <- data.frame('wl'= wl,'fwhm'=fwhm)                             #wl group by fwhm in data.frame
  SRF <- prosail::GetRadiometry(SpectralProps = SpectralProps,
                                SensorName = paste('Filter_',fwhm,sep=''),
                                Path_SensorResponse = SRF_Dir,
                                SaveSRF = FALSE)
  # Inversion
  Parms2Estimate <- "ALL"
  SpecPROSPECT_Sensor <- data.frame(prosail::applySensorCharacteristics(wvl = SpecPROSPECT$lambda,
                                                                        InRefl = SpecPROSPECT,
                                                                        SRF = SRF))
  RT_filter <- gaussian_filter(SRF = SRF, lambda = lambda, Refl=Refl, Tran=Tran)
  nbSamples <- ncol(Refl)
  message(paste('FWHM = ',fwhm))
  parms_est_df <- data.frame()
  pb <- progress::progress_bar$new(
    format = "PROSPECT inversion [:bar] :percent in :elapsedfull",
    total = nbSamples, clear = FALSE, width= 100)
  for (i in 1:nbSamples) {                                         #filtered reflectance-transmittance list
    pb$tick()
    Rtmp <- RT_filter$Refl[,i]
    Ttmp <- RT_filter$Tran[,i]
    Invert_est <- Invert_PROSPECT(SpecPROSPECT =  SpecPROSPECT_Sensor,
                                  Refl = Rtmp,
                                  Tran = Ttmp,
                                  Parms2Estimate = Parms2Estimate,
                                  PROSPECT_version = 'D')
    parms_est_df <- rbind(parms_est_df, Invert_est)
  }
  return(parms_est_df)
}


