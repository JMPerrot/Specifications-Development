################################################################################
## This code makes it possible to test different step values for the prospect
## model in order to follow the evolution of the precision of the model step by
## step. There is no filter
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

################################################################################
# input output directories
################################################################################
PathData <- '../../01_DATA'
PathResults <- '../../03_RESULTS'
PAS_Dir <- file.path(PathResults,'01_TEST_PAS_MAT')
dir.create(path = PAS_Dir,showWarnings = F,recursive = T)
################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)

################################################################################
# output directories for fig
################################################################################

for (i in c("LMA","EWT","N")){
  PAS_Dir_Fig <- file.path(PAS_Dir,'FIGURES',i)
  dir.create(path = PAS_Dir_Fig,showWarnings = F,recursive = T)
  
  
  List_pas<-c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,175,200,225,250)
  
  for (pas in List_pas) {
    if (!file.exists(paste(PAS_Dir,paste("/PAS=",pas,"nm.RData", sep=""), sep = ""))){
      ############################################################################
      # definition of the sampling step
      ############################################################################
      wl<-seq(1000,2450, by = pas)
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
      T1<-Sys.time()
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
      
      T2<-Sys.time()
      Time<-as.numeric(T2-T1)
      TIME<-data.frame("TIME"=Time)
      parms_est_df<-cbind(parms_est_df,TIME)
      
      filename<-paste(PAS_Dir,paste("/PAS=",pas,"nm.RData", sep=""), sep = "")
      save(parms_est_df,file = filename)
    }else{
      load(paste(PAS_Dir,paste("/PAS=",pas,"nm.RData", sep=""), sep = ""))
    }
    
    fileName = paste(PAS_Dir_Fig,paste("/PAS=",pas,"nm_",i,".png", sep=""), sep = "")
    
    if(i=="LMA"){
      scatter_inversion(target = Biochemistry$LMA,
                        estimate = parms_est_df$LMA,
                        Colors = "red", 
                        fileName = fileName, 
                        Labs = list("Mesured LMA","Estimated LMA"),
                        PlotStats = TRUE,
                        categories = "RT")
    }else if(i=="EWT"){
      scatter_inversion(target = Biochemistry$EWT,
                        estimate = parms_est_df$EWT,
                        Colors = "red", 
                        fileName = fileName, 
                        Labs = list("Mesured EWT","Estimated EWT"),
                        PlotStats = TRUE,
                        categories = "RT")
    }
  }
}

