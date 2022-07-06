################################################################################
## This code aims at exploring the influence of spectral sampling on the 
## performances of PROSPECT inversion
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
library(doFuture)
source('../../Libraries/Lib_Plots.R')
source('../../Libraries/Lib_Analysis_Inversion.R')

################################################################################
# input output directories
################################################################################
PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/R&T'
SpectralSampling_Dir <- file.path(PathResults,'02_SpectralSampling')
dir.create(path = SpectralSampling_Dir,
           showWarnings = F,recursive = T)

################################################################################
# load leaf optics dataset
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
lambda <- Reflectance$V1
Reflectance <- Reflectance[,-1]
Transmittance <- Transmittance[,-1]

################################################################################
# Define parameters for inversion
################################################################################
Parms2Estimate <- c('EWT','LMA','CHL','CAR')
Parms2Estimate_ind <- list('EWT'=c('EWT'),"LMA"=c('LMA'),"CHL"=c('CHL'),"CAR"=c('CAR'))
# define spectral sampling
SpectralSampling <- list()
SpectralSampling$CHL<-SpectralSampling$CAR <- as.list(c(seq(1,50,by=1)))
SpectralSampling$EWT<-SpectralSampling$LMA <- as.list(c(seq(2,50,by=1)))

# define spectral domain
SpectralDomain <- list()
SpectralDomain$CHL<-SpectralDomain$CAR <- list('minWL' = 400, 'maxWL' = 900)
SpectralDomain$EWT<-SpectralDomain$LMA <- list('minWL' = 1300, 'maxWL' = 2400)

# define parameters to estimate during inversion
ParmsEstInv <- list()
ParmsEstInv$CAR <- ParmsEstInv$CHL <- "ALL"
ParmsEstInv$EWT <- ParmsEstInv$LMA <- c('EWT', 'LMA', 'N')

# define parameters to estimate during inversion
InitValues <- list()
InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)


################################################################################
# Perform inversion
################################################################################
# for each parameter
nbWorkers <- 4
registerDoFuture()
plan(multisession, workers = nbWorkers)
Estimate_SpectralSampling <- list()
SpectralSampling_subDir <- file.path(SpectralSampling_Dir,'_individualSteps')
dir.create(path = SpectralSampling_subDir,
           showWarnings = F,recursive = T)

for (parm in Parms2Estimate){
  print(parm)
  # apply multiprocessing for each spectral sampling
  Invert_SpectralSampling <- function() {
    foreach(subsample = SpectralSampling[[parm]]) %dopar% {
      Invert_est<-list()
      if (subsample<10){
        subChar <- paste('00',as.character(subsample),sep = '')
      }else if (subsample<100){
        subChar <- paste('0',as.character(subsample),sep = '')
      }else if (subsample<1000){
        subChar <- as.character(subsample)
      }
      FileName <- file.path(SpectralSampling_subDir,
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
        
        Invert_est <- prospect::Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                                Refl = DataFit$Refl,
                                                Tran = DataFit$Tran,
                                                PROSPECT_version = 'D',
                                                Parms2Estimate = ParmsEstInv[[parm]],
                                                InitValues = InitValues)
        save(Invert_est,file = FileName)
      } else {
        load(file = FileName)
      }
      return(Invert_est)
    }
  }
  PROSPECT_EST <- Invert_SpectralSampling()
  for (subparm in Parms2Estimate_ind[[parm]]){
    Estimate_SpectralSampling[[subparm]] <- do.call(cbind,lapply(PROSPECT_EST,'[[',subparm))
    Estimate_SpectralSampling[[subparm]] <- data.frame(Estimate_SpectralSampling[[subparm]])
    colnames(Estimate_SpectralSampling[[subparm]]) <- unlist(SpectralSampling[[parm]])
  }
}
plan(sequential)

# save all spectral sampling results in a unique file per parameter
for (parm in names(Estimate_SpectralSampling)){
  FileName <- file.path(SpectralSampling_Dir,paste(parm,'_SpecSampling.csv',sep = ''))
  write_delim(x = Estimate_SpectralSampling[[parm]],
              file = FileName,
              delim = '\t',
              col_names = T)
}

# ################################################################################
# # Produce figures
# ################################################################################

# plot results
PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
MinMax <- list('CHL' = c(0,120), 'CAR'  = c(0,30), 'LMA' = c(0,40), 'EWT' = c(0,60))
UnitsParms <- list('CHL'='(µg/cm²)', 'CAR'='(µg/cm²)', 'EWT'='(mg/cm²)', 'LMA'='(mg/cm²)')
Factor <- list('CHL'=1, 'CAR'=1, 'EWT'=1000, 'LMA'=1000)
PlotObj <- list()
for (parm in names(Estimate_SpectralSampling)){
  Labs <- c(paste('Measured',parm,UnitsParms[[parm]]), 
            paste('Estimated',parm,UnitsParms[[parm]]))
  
  for (col0 in colnames(Estimate_SpectralSampling[[parm]])){
    print(c(parm,col0))
    valStep <- as.numeric(col0)
    if (valStep<10){subChar <- paste('00',as.character(valStep),sep = '')}
    else if (valStep<100){subChar <- paste('0',as.character(valStep),sep = '')}
    else if (valStep<1000){subChar <- as.character(valStep)}
    fileName <- file.path(SpectralSampling_subDir,'FIGURES', paste(parm,'_',subChar,'.png',sep = ''))
    PlotObj <- scatter_inversion(target = Factor[[parm]]*Biochemistry[[parm]],
                                 estimate = Factor[[parm]]*Estimate_SpectralSampling[[parm]][[col0]],
                                 Colors = PlotCols[[parm]],
                                 Labs = Labs,
                                 fileName = fileName,
                                 categories = "RT_FS",
                                 PlotStats = T,
                                 MinMaxAxis = MinMax[[parm]],
                                 size = 1)
  }
}
