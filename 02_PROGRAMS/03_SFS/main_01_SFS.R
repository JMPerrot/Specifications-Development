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
PAS_Dir <- file.path(PathResults,'02_TEST_PAS')
dir.create(path = PAS_Dir,showWarnings = F,recursive = T)
################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
################################################################################

parms_est_df=data.frame("CAR"=NA, "deleted_wl" = NA , "obs" = NA)              #           

seqWL<-c(500:540)#seq(400, 1000, by =10)


while(length(seqWL)>5){
  nrmse_df = data.frame("deleted_wl"=NA, "NRMSE" = NA)
  parms_est_df=data.frame("CAR"=NA, "deleted_wl" = NA , "obs" = NA) 
  for (i in c(1:length(seqWL))){
    seqWLR<-seqWL[-i] #sequence de longueure d'onde reduite
    lambda<-c(unlist(Reflectance[,1]))
    WLselect<-match(seqWLR,lambda)
    
    pb <- progress::progress_bar$new(
      format = "PROSPECT inversion [:bar] :percent in :elapsedfull",
      total = length(Reflectance), clear = FALSE, width= 100)
    
    for (j in c(2:length(Reflectance))) {
      DataFit<-FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                               lambda = SpecPROSPECT$lambda,
                               Refl = matrix(unlist(Reflectance[[j]])),
                               Tran = matrix(unlist(Transmittance[[j]])),
                               UserDomain = seqWLR,
                               UL_Bounds = FALSE)
      
      Invert_est<-Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                  Refl = DataFit$Refl,
                                  Tran = DataFit$Tran,
                                  InitValues = data.frame(CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0.01, EWT = 0.01, LMA
                                                          = 0.01, PROT = 0.001, CBC = 0.009, N = 1.5, alpha = 40),
                                  Parms2Estimate = c('CHL','CAR'))
      
      est = data.frame("CAR"=Invert_est$CAR, 
                       "deleted_wl" = seqWL[i],
                       "obs" = j)    
      
      parms_est_df=rbind(parms_est_df, est)
      pb$tick()
    }
    parms<- parms_est_df %>% filter(parms_est_df$deleted_wl == seqWL[i])
    nrmse_wl=get_performances_inversion(target = Biochemistry$CAR,
                                        estimate = parms$CAR,
                                        categories= TRUE)[3]
    nrmse_df=rbind(nrmse_df, data.frame("deleted_wl" = seqWL[i],"NRMSE" = nrmse_wl[[1]]))
  }
  
  nrmse_df<-nrmse_df[-1,]
  min_nrmse<-min(as.matrix(nrmse_df$NRMSE))
  nline<-match(1,match(nrmse_df$NRMSE,min_nrmse))
  seqWL<-seqWL[-nline]
}

nrmse_df<-nrmse_df[-1,]
