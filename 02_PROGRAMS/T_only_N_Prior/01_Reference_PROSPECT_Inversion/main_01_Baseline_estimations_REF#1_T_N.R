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
library(tidyverse)
library(prospect)
library(prosail)
library(data.table)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
source('../../Libraries/Lib_Analysis_Inversion.R')
source('../../Libraries/Lib_Plots.R')

################################################################################
# input output directories
################################################################################
PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/T_only_N_Prior/01_Reference'
dir.create(PathResults,showWarnings = F,recursive = T)
################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
# Get the wavelengths corresponding to the reflectance and transmittance measurements  
lambda <- Reflectance[[1]]
Refl<- Reflectance[,-1]
Tran <- Transmittance[,-1]

# Get the number of samples
nbSamples <- ncol(Refl)

# Estimate all parameters for PROSPECT-D
Parms2Estimate<-list()
Parms2Estimate$CHL<-Parms2Estimate$CAR  <- "ALL"
Parms2Estimate$LMA<-Parms2Estimate$EWT  <- c('EWT','LMA')

N_prior <- Get_Nprior(SpecPROSPECT = SpecPROSPECT, 
                      lambda = lambda, 
                      Refl = NULL, 
                      Tran = Tran,
                      OptWL_R = list(NIR = 800, SWIR = 1131),  
                      OptWL_T = list(NIR = 753, SWIR = 1121))
InitValues<-list()
for (i in c(1:nbSamples)) {
  InitValues[[i]] <- data.frame("CHL"=40, "CAR"=10, "ANT"=0.1, "BROWN"=0, "EWT"=0.01, "LMA"=0.01, "N"=N_prior[i,])
}
# Adjust spectral domain for SpecPROSPECT to fit leaf optical properties 
SubData<-list()
SubSpecPROSPECT<-list()
Sublambda<-list()
SubRefl<-list()
SubTran<-statsInv<- res <- res1 <- list()

Inversion_Ref1 <- R2_Refl1 <- NRMSE_Refl1 <- list()
SubData$CHL<-SubData$CAR <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,
                                            lambda = lambda, Refl = NULL, Tran =Tran,
                                            UserDomain = c(400,900),
                                            UL_Bounds = TRUE)
SubData$LMA<-SubData$EWT <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,
                                            lambda = lambda, Refl = NULL, Tran =Tran,
                                            UserDomain = c(1300,2400),
                                            UL_Bounds = TRUE)
for (parm in c("CHL","CAR","LMA","EWT")) {
  SubSpecPROSPECT[[parm]] <- SubData[[parm]]$SpecPROSPECT
  Sublambda[[parm]] <- SubData[[parm]]$lambda
  SubRefl[[parm]] <- SubData[[parm]]$Refl
  SubTran[[parm]] <- SubData[[parm]]$Tran
  
  print('PROSPECT inversion using full spectral range')
  for (i in c(1:nbSamples)) {
    res[[parm]][[i]] <- Invert_PROSPECT(SpecPROSPECT = SubSpecPROSPECT[[parm]],
                                   Refl = NULL, 
                                   Tran = SubTran[[parm]][[i]],
                                   PROSPECT_version = 'D', 
                                   Parms2Estimate = Parms2Estimate[[parm]],
                                   InitValues = InitValues[[i]])
  }
  for (i in c(1:nbSamples)) {
    res1[[parm]]<-rbind(res1[[parm]],res[[parm]][[i]][[parm]])
  }
  # compute statistics for inversion
  ParmsOfInterest <- c('CHL', 'CAR', 'EWT', 'LMA')
  UnitsParms <- list('CHL'='(µg/cm²)', 'CAR'='(µg/cm²)', 'EWT'='(mg/cm²)', 'LMA'='(mg/cm²)')
  Factor <- list('CHL'=1, 'CAR'=1, 'EWT'=1000, 'LMA'=1000)
  
  
  Inversion_Ref1[[parm]][[parm]] <- data.frame('measured' = Factor[[parm]]*Biochemistry[[parm]],
                                               'estimated' = Factor[[parm]]*res1[[parm]])#res[[parm]][[parm]])
  
  statsInv[[parm]] <- get_performances_inversion(target = Inversion_Ref1[[parm]][[parm]][["measured"]],
                                                 estimate = Inversion_Ref1[[parm]][[parm]][["estimated"]],
                                                 categories = T)
  
  R2_Refl1[[parm]] <- statsInv$R2
  NRMSE_Refl1[[parm]] <- statsInv$NRMSE
  
}
# save results
for (parm in c("CHL","CAR","LMA","EWT")){
  FileName <- file.path(PathResults,paste(parm,'_REFERENCE#1_T_N.RData',sep = ''))
  ResultsInversion <- Inversion_Ref1[[parm]][[parm]] 
  save(ResultsInversion ,file = FileName)
}

# plot results
PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
MinMax <- list('CHL' = c(0,120), 'CAR'  = c(0,30), 'LMA' = c(0,40), 'EWT' = c(0,60))
PlotObj <- list()
for (parm in c("CHL","CAR","LMA","EWT")){
  Labs <- c(paste('Measured',parm,UnitsParms[[parm]]), 
            paste('Estimated',parm,UnitsParms[[parm]]))
  fileName <- file.path(PathResults,'FIGURES', paste(parm,'_REFERENCE#1_T_N.png',sep = ''))
  PlotObj[[parm]] <- scatter_inversion(target = Inversion_Ref1[[parm]][[parm]]$measured,
                                       estimate = Inversion_Ref1[[parm]][[parm]]$estimated,
                                       Colors = PlotCols[[parm]],
                                       Labs = Labs,
                                       fileName = fileName,
                                       categories = "T_N",
                                       PlotStats = T, 
                                       MinMaxAxis = MinMax[[parm]],
                                       size = 1)
}
gg <- grid.arrange(grobs = PlotObj, ncol = 2)
fileName <- file.path(PathResults,'FIGURES', 'All_REFERENCE#1_T_N.png')
ggsave(fileName, gg,device = "png")
