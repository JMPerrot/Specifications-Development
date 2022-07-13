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
library(doFuture)
source('../../Libraries/Lib_Analysis_Inversion.R')
source('../../Libraries/Lib_Plots.R')

################################################################################
# input output directories
################################################################################
PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/T_only_N_Prior/02_SpecSampling'
dir.create(PathResults,showWarnings = F,recursive = T)

################################################################################
# load leaf optics dataset
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
lambda <- Reflectance$V1
Reflectance <- Reflectance[,-1]
Transmittance<- Transmittance[,-1]
nbSamples <- ncol(Reflectance)

################################################################################
# Define parameters for inversion
################################################################################
Parms2Estimate <- c('EWT_LMA','CHL_CAR')
Parms2Estimate_ind <- list('EWT_LMA'=c('EWT','LMA'),'CHL_CAR'=c('CHL','CAR'))
# define spectral sampling
SpecSampling <- list()
SpecSampling$CHL_CAR <- as.list(c(seq(2,50,by=1)))
SpecSampling$EWT_LMA <- as.list(c(seq(2,50,by=1)))

# define spectral domain
SpectralDomain <- list()
SpectralDomain$CHL_CAR <- list('minWL' = 400, 'maxWL' = 900)
SpectralDomain$EWT_LMA <- list('minWL' = 1300, 'maxWL' = 2400)

# define parameters to estimate during inversion
ParmsEstInv <- list()
ParmsEstInv$CHL_CAR <- c('CHL', 'CAR','ANT', 'EWT', 'LMA')
ParmsEstInv$EWT_LMA <- c('EWT', 'LMA')
N_prior <- Get_Nprior(SpecPROSPECT = SpecPROSPECT, lambda = lambda, Refl = Reflectance, Tran = NULL)


# define parameters to estimate during inversion
InitValues <- list()
InitValues<-list()
for (i in c(1:nbSamples)) {
  InitValues[[i]] <- data.frame("CHL"=40, "CAR"=10, "ANT"=0.1, "BROWN"=0, "EWT"=0.01, "LMA"=0.01, "N"=N_prior[i,])
}

################################################################################
# Perform inversion
################################################################################
# for each parameter
nbWorkers <- 2
registerDoFuture()
plan(multisession, workers = nbWorkers)
Estimate_SpecSampling <- list()
SpecSampling_subDir <- file.path(PathResults,'_individualSteps')
dir.create(path = SpecSampling_subDir,
           showWarnings = F,recursive = T)

for (parm in Parms2Estimate){
  print(parm)
  # apply multiprocessing for each spectral sampling
  Invert_SpecSampling <- function() {
    foreach(subsample = SpecSampling[[parm]]) %dopar% {
      if (subsample<10){
        subChar <- paste('00',as.character(subsample),sep = '')
      }else if (subsample<100){
        subChar <- paste('0',as.character(subsample),sep = '')
      }else if (subsample<1000){
          subChar <- as.character(subsample)
          }
      FileName <- file.path(SpecSampling_subDir,
                            paste(parm,'_SpecSampling_',subChar,'.RData',sep = ''))
      if (!file.exists(FileName)){
        lambda_tmp <- seq(from = SpectralDomain[[parm]]$minWL,
                          to = SpectralDomain[[parm]]$maxWL,
                          by = subsample)
        
        DataFit <- FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                                   lambda = lambda,
                                   Refl = Reflectance,
                                   Tran = Transmittance,
                                   UserDomain = lambda_tmp,
                                   UL_Bounds = FALSE)
        Invert_est<-list()
        for (i in c(1:nbSamples)) {
          Invert <- prospect::Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                                  Refl = NULL,
                                                  Tran = DataFit$Tran[[i]],
                                                  PROSPECT_version = 'D',
                                                  Parms2Estimate = ParmsEstInv[[parm]],
                                                  InitValues = InitValues[[i]])
          Invert_est<-rbind(Invert_est,Invert)
        }
        save(Invert_est,file = FileName)
      } else {
        load(file = FileName)
      }
      return(Invert_est)
    }
  }
  PROSPECT_EST <- Invert_SpecSampling()
  for (subparm in Parms2Estimate_ind[[parm]]){
    Estimate_SpecSampling[[subparm]] <- do.call(cbind,lapply(PROSPECT_EST,'[[',subparm))
    Estimate_SpecSampling[[subparm]] <- data.frame(Estimate_SpecSampling[[subparm]])
    colnames(Estimate_SpecSampling[[subparm]]) <- unlist(SpecSampling[[parm]])
  }
}
plan(sequential)

# save all spectral sampling results in a unique file per parameter
for (parm in names(Estimate_SpecSampling)){
  FileName <- file.path(PathResults,paste(parm,'_SpecSampling.csv',sep = ''))
  write_delim(x = Estimate_SpecSampling[[parm]],
              file = FileName,
              delim = '\t',
              col_names = T)
}

# ################################################################################
# # Produce figures
# ################################################################################

  # plot results
PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
MinMax <- list('CHL' = c(0,120), 'CAR'  = c(0,30), 'LMA' = c(0,65), 'EWT' = c(0,105))
UnitsParms <- list('CHL'='(µg/cm²)', 'CAR'='(µg/cm²)', 'EWT'='(mg/cm²)', 'LMA'='(mg/cm²)')
Factor <- list('CHL'=1, 'CAR'=1, 'EWT'=1000, 'LMA'=1000)
PlotObj <- list()
for (parm in names(Estimate_SpecSampling)){
  Labs <- c(paste('Measured',parm,UnitsParms[[parm]]), 
            paste('Estimated_T_Only_N_priori',parm,UnitsParms[[parm]]))
  
  for (col0 in colnames(Estimate_SpecSampling[[parm]])){
    print(c(parm,col0))
    valStep <- as.numeric(col0)
    if (valStep<10){subChar <- paste('00',as.character(valStep),sep = '')}
    else if (valStep<100){subChar <- paste('0',as.character(valStep),sep = '')}
    else if (valStep<1000){subChar <- as.character(valStep)}
    fileName <- file.path(SpecSampling_subDir,'FIGURES', paste(parm,'_',subChar,'.png',sep = ''))
    PlotObj <- scatter_inversion(target = Factor[[parm]]*Biochemistry[[parm]],
                                 estimate = Factor[[parm]]*Estimate_SpecSampling[[parm]][[col0]],
                                 Colors = PlotCols[[parm]],
                                 Labs = Labs,
                                 fileName = fileName,
                                 categories = "T_N",
                                 PlotStats = T,
                                 MinMaxAxis = MinMax[[parm]],
                                 size = 1)
  }
}



