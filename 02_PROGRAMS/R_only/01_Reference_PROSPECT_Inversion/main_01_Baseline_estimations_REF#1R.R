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
library(data.table)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
source('../Libraries/Lib_Analysis_Inversion.R')
source('../Libraries/Lib_Plots.R')

################################################################################
# input output directories
################################################################################
PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/01_Reference/R_only'
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
Parms2Estimate  = 'ALL'
InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
# Adjust spectral domain for SpecPROSPECT to fit leaf optical properties 
SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,
                           lambda = lambda, Refl = Refl, Tran =Tran,
                           UserDomain = c(lambda[1],tail(lambda,n = 1)),
                           UL_Bounds = TRUE)
SubSpecPROSPECT <- SubData$SpecPROSPECT
Sublambda <- SubData$lambda
SubRefl <- SubData$Refl
SubTran <- SubData$Tran

print('PROSPECT inversion using full spectral range')
res <- Invert_PROSPECT(SpecPROSPECT = SubSpecPROSPECT, 
                       Refl = SubRefl, 
                       Tran = SubTran, 
                       PROSPECT_version = 'D', 
                       Parms2Estimate = Parms2Estimate, 
                       InitValues = InitValues)

# compute statistics for inversion
ParmsOfInterest <- c('CHL', 'CAR', 'EWT', 'LMA')
UnitsParms <- list('CHL'='(µg/cm²)', 'CAR'='(µg/cm²)', 'EWT'='(mg/cm²)', 'LMA'='(mg/cm²)')
Factor <- list('CHL'=1, 'CAR'=1, 'EWT'=1000, 'LMA'=1000)
Inversion_Ref1 <- R2_Refl1 <- NRMSE_Refl1 <- list()
for (parm in ParmsOfInterest){
  Inversion_Ref1[[parm]] <- data.frame('measured' = Factor[[parm]]*Biochemistry[[parm]],
                                       'estimated' = Factor[[parm]]*res[[parm]])
  
  statsInv <- get_performances_inversion(target = Inversion_Ref1[[parm]]$measured,
                                         estimate = Inversion_Ref1[[parm]]$estimated,
                                         categories = T)
  
  R2_Refl1[[parm]] <- statsInv$R2
  NRMSE_Refl1[[parm]] <- statsInv$NRMSE
}

# save results
for (parm in ParmsOfInterest){
  FileName <- file.path(PathResults,paste(parm,'_REFERENCE#1.RData',sep = ''))
  ResultsInversion <- Inversion_Ref1[[parm]]
  save(ResultsInversion ,file = FileName)
}
  
# plot results
PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
MinMax <- list('CHL' = c(0,120), 'CAR'  = c(0,30), 'LMA' = c(0,40), 'EWT' = c(0,60))
PlotObj <- list()
for (parm in ParmsOfInterest){
  Labs <- c(paste('Measured',parm,UnitsParms[[parm]]), 
            paste('Estimated',parm,UnitsParms[[parm]]))
  fileName <- file.path(PathResults,'FIGURES', paste(parm,'_REFERENCE#1.png',sep = ''))
  PlotObj[[parm]] <- scatter_inversion(target = Inversion_Ref1[[parm]]$measured,
                                       estimate = Inversion_Ref1[[parm]]$estimated,
                                       Colors = PlotCols[[parm]],
                                       Labs = Labs,
                                       fileName = fileName,
                                       categories = "RT_FS",
                                       PlotStats = T, 
                                       MinMaxAxis = MinMax[[parm]],
                                       size = 1)
}
gg <- grid.arrange(grobs = PlotObj, ncol = 2)
fileName <- file.path(PathResults,'FIGURES', 'All_REFERENCE#1.png')
ggsave(fileName, gg,device = "png")
