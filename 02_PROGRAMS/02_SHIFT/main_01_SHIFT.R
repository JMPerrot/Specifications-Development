################################################################################
## This code makes it possible to test different offset/shift values for the 
## prospect model in order to follow the evolution of the precision of the 
## model. There is no filter
################################################################################

# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
################################################################################
# Libraries required
################################################################################
library(tidyverse)
library(prospect)
library(data.table)
source('../Libraries/Lib_Plots.R')
source('../Libraries/Lib_Analysis_Inversion.R')
source('../Libraries/Lib_Inversion_PROSPECT.R')

################################################################################
# input output directories
################################################################################
PathData <- '../../01_DATA'
PathResults <- '../../03_RESULTS'
SHIFT_Dir <- file.path(PathResults,'02_TEST_SHIFT')
dir.create(path = SHIFT_Dir,showWarnings = F,recursive = T)

################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)

pas <- 100
shift_list <- c(0:pas)

SHIFT_Dir_Fig <- file.path(SHIFT_Dir,paste('FIGURES/CAR/CAR_',pas,'nm', sep =""))
dir.create(path = SHIFT_Dir_Fig,showWarnings = F,recursive = T)
#SHIFT_Dir_Data<-paste(SHIFT_Dir,paste("/pas=",pas,"nm/SHIFT=",shift,"nm",".RData", sep=""), sep = "")
SHIFT_Dir_Data<-paste(SHIFT_Dir,paste("/pas=",pas,"nm", sep=""), sep = "")
dir.create(paste(SHIFT_Dir,paste("/pas=",pas,"nm", sep=""), sep = ""),showWarnings = F,recursive = T)

for (shift in shift_list) {
  if (!file.exists(paste(SHIFT_Dir,paste("/pas=",pas,"nm/SHIFT=",shift,"nm",".RData", sep=""), sep = ""))){
    ############################################################################
    # definition of the sampling step
    ############################################################################
    wl<-seq(400+shift,1000, by = pas)
    lamb <- c(unlist(Reflectance[,1]))
    WLselect<-match(wl,lamb)
    lambda<-lamb[WLselect]
    ############################################################################
    # Estimate all parameters for PROSPECT-D
    ############################################################################
    parms_est_df=data.frame("CHL" = NA, "CAR" = NA, "ANT" = NA,
                            "EWT" = NA, "LMA" = NA, "N" = NA)
    Parms2Estimate  = c('CHL','CAR','ANT','EWT','LMA') 
    #we are estimating all parameters but spectral band is not optimized for LMA and EWT
    
    InitValues <- data.frame(CHL=45, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
    
    pb <- progress::progress_bar$new(
      format = "PROSPECT inversion [:bar] :percent in :elapsedfull",
      total = length(Reflectance), clear = FALSE, width= 100)
    
    # parms_est_df <- lapply( X = c(2:length(Reflectance)), FUN =  Invert_PROSPECT_USERDOMAIN,
    #                        SpecPROSPECT = SpecPROSPECT,
    #                        Reflectance = Reflectance, 
    #                        Transmittance= Transmittance,
    #                        lbd = lambda,
    #                        parms_est_df = parms_est_df)
    for (j in c(1:length(Reflectance))) {
    DataFit<-FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                             lambda = SpecPROSPECT$lambda,
                             Refl = matrix(unlist(Reflectance[[j]])),
                             Tran = matrix(unlist(Transmittance[[j]])),
                             UserDomain = lambda,
                             UL_Bounds = FALSE)

    Invert_est <- prospect::Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                            Refl = DataFit$Refl,
                                            Tran = DataFit$Tran,
                                            PROSPECT_version = 'D',
                                            Parms2Estimate = Parms2Estimate)
    est=data.frame("CHL" = Invert_est$CHL,"CAR" = Invert_est$CAR,"ANT" = Invert_est$ANT,
                   "EWT" = Invert_est$EWT, "LMA" = Invert_est$LMA, "N" = Invert_est$N)
    parms_est_df=rbind(parms_est_df, est)
    pb$tick()
    }
    parms_est_df<-parms_est_df[-c(1,2),]
    filename<-paste(SHIFT_Dir,paste("/pas=",pas,"nm/SHIFT=",shift,"nm",".RData", sep=""), sep = "")
    save(parms_est_df,file = filename)
  }else{
    load(paste(SHIFT_Dir,paste("/pas=",pas,"nm/SHIFT=",shift,"nm",".RData", sep=""), sep = ""))
  }
  
  fileName = paste(SHIFT_Dir_Fig,paste("/SHIFT=",shift,"nm_CAR.png", sep=""), sep = "")
  scatter_inversion(target = Biochemistry$CAR,
                    estimate = parms_est_df$CAR,
                    Colors = "red", 
                    fileName = fileName, 
                    Labs = list("Mesured CAR","Estimated CAR"),
                    PlotStats = TRUE,
                    categories = "RT")
}
