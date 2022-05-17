# ==============================================================================
# Librairie test pour applique Lapply 
# NE FONCTIONNE PAS POUR LE MOMENT
# ==============================================================================
# PROGRAMMERS: Jean-Baptiste FERET, Jean-Malo PERROT
#
# ==============================================================================
# This Library includes nmrse_func
# ==============================================================================

#' This function allows application of a gaussian filter on reflectance and
#' transmittance data
#'
#' @param wl int : .
#'
#'
#' @return nrmse : Normalized RMSE
#' @importFrom 
#' @export
#' @importFrom
#' @export
#' 
Invert_PROSPECT_USERDOMAIN <-function(SpecPROSPECT, Reflectance, Transmittance,wl,lbd, parms_est_df){
    DataFit<-FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                             lambda = SpecPROSPECT$lambda,
                             Refl = matrix(unlist(Reflectance[[wl]])),
                             Tran = matrix(unlist(Transmittance[[wl]])),
                             UserDomain = lbd,
                             UL_Bounds = FALSE)
    
    Invert_est <- prospect::Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                            Refl = DataFit$Refl,
                                            Tran = DataFit$Tran,
                                            PROSPECT_version = 'D',
                                            Parms2Estimate = c('CHL','CAR','ANT','EWT','LMA'))
    
    est=data.frame("CHL" = Invert_est$CHL,"CAR" = Invert_est$CAR,"ANT" = Invert_est$ANT,
                   "EWT" = Invert_est$EWT, "LMA" = Invert_est$LMA, "N" = Invert_est$N)
    parms_est_df=rbind(parms_est_df, est)
    # parms_est_df<-parms_est_df[-c(1,2),]
    filename<-paste(SHIFT_Dir,paste("/pas=",pas,"nm/SHIFT=",shift,"nm",".RData", sep=""), sep = "")
    save(parms_est_df,file = filename)
    pb$tick()
    return(parms_est_df)
}