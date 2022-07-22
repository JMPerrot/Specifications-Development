################################################################################
## This code aims at exploring the influence of spectral sampling on the 
## performances of PROSPECT inversion
################################################################################
# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Libraries required -----------------------------------------------------------

library(tidyverse)
library(prospect)
library(data.table)
library(doFuture)
source('../../Libraries/Lib_Plots.R')
source('../../Libraries/Lib_Analysis_Inversion.R')


# input output directories -----------------------------------------------------

PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/R&T'
SpectralShifting_Dir <- file.path(PathResults,'03_SpectralShifting')
SpectralSampling_Dir <- file.path(PathResults,'02_SpectralSampling')
dir.create(path = SpectralShifting_Dir, showWarnings = F,recursive = T)






# load leaf optics dataset -----------------------------------------------------

dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
lambda <- Reflectance$V1
Reflectance <- Reflectance[,-1]
Transmittance <- Transmittance[,-1]


# Define parameters for inversion ----------------------------------------------

Parms2Estimate <- c('EWT','LMA','CHL','CAR')
Parms2Estimate_ind <- list('CHL'=c('CHL'),'EWT'= 'EWT','LMA'=c('LMA'),'CAR'=c('CAR'))

Stats<-list()
Opt_Sampling<-list()
for (parm in Parms2Estimate){
  FileName <- file.path(SpectralSampling_Dir,paste(parm,'_Statistics.csv',sep = ''))
  Stats[[parm]] <- readr::read_delim(file = FileName,delim = '\t')
  Opt_Sampling[[parm]] <- Stats[[parm]]$Sampling[which(Stats[[parm]]$NRMSE==min(Stats[[parm]]$NRMSE))]
}
# define spectral sampling
SpectralSampling <- list()
SpectralSampling$CHL <- Opt_Sampling$CHL 
SpectralSampling$CAR <- Opt_Sampling$CAR
SpectralSampling$LMA <- 45
SpectralSampling$EWT <- 45
# define spectral domain
minlambda <- list()
minlambda$CHL <- c(0:max(SpectralSampling[["CHL"]]))
minlambda$EWT <- c(0:max(SpectralSampling[["EWT"]]))
minlambda$CAR <- c(0:max(SpectralSampling[["CAR"]]))
minlambda$LMA <- c(0:max(SpectralSampling[["LMA"]]))


# define parameters to estimate during inversion
ParmsEstInv <- list()
ParmsEstInv$CHL<-ParmsEstInv$CAR <- "ALL"
ParmsEstInv$EWT<-ParmsEstInv$LMA <- c('EWT', 'LMA', 'N')

# define parameters to estimate during inversion
InitValues <- list()
InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)


# Perform inversion ------------------------------------------------------------

# for each parameter
nbWorkers <- 4 # number of core
registerDoFuture()
plan(multisession, workers = nbWorkers)
Estimate_SpectralShifting <- list()
SpectralShifting_subDir <- file.path(SpectralShifting_Dir,'_individualSteps')
dir.create(path = SpectralShifting_subDir,
           showWarnings = F,recursive = T)

for (parm in Parms2Estimate){
  print(parm)
  # apply multiprocessing for each spectral sampling
  Invert_SpectralShifting <- function() {
    foreach(subshifting = c(0:max(unlist(minlambda)))) %dopart% {
      SpectralDomain <- Invert_est <- list()
      SpectralDomain$CHL<-SpectralDomain$CAR <- list('minWL' = 400 + subshifting, 'maxWL' = 900) 
      SpectralDomain$EWT<-SpectralDomain$LMA <- list('minWL' = 1300+subshifting, 'maxWL' = 2400)
      
      if (subshifting<10){
        subChar <- paste('00',as.character(subshifting),sep = '')
      }else if (subshifting<100){
        subChar <- paste('0',as.character(subshifting),sep = '')
      }else if (subshifting<1000){
          subChar <- as.character(subshifting)
          }
      FileName <- file.path(SpectralShifting_subDir,
                            paste(parm,'_SpecShifting_',subChar,'.RData',sep = ''))
      if (!file.exists(FileName)){
        lambda_tmp <- seq(from = SpectralDomain[[parm]]$minWL,
                          to = SpectralDomain[[parm]]$maxWL,
                          by = SpectralSampling[[parm]])
        
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
  PROSPECT_EST <- Invert_SpectralShifting()
  for (subparm in Parms2Estimate_ind[[parm]]){
    Estimate_SpectralShifting[[subparm]] <- do.call(cbind,lapply(PROSPECT_EST,'[[',subparm))
    Estimate_SpectralShifting[[subparm]] <- data.frame(Estimate_SpectralShifting[[subparm]])
    colnames(Estimate_SpectralShifting[[subparm]]) <- minlambda[[parm]]
  }
}
plan(sequential)

# save all spectral sampling results in a unique file per parameter
for (parm in names(Estimate_SpectralShifting)){
  FileName <- file.path(SpectralShifting_Dir,paste(parm,'_SpecShifting.csv',sep = ''))
  write_delim(x = Estimate_SpectralShifting[[parm]],
              file = FileName,
              delim = '\t',
              col_names = T)
}


## Produce figures ------------------------------------------------------------


# plot results
PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
MinMax <- list('CHL' = c(0,120), 'CAR'  = c(0,30), 'LMA' = c(0,40), 'EWT' = c(0,60))
UnitsParms <- list('CHL'='(µg/cm²)', 'CAR'='(µg/cm²)', 'EWT'='(mg/cm²)', 'LMA'='(mg/cm²)')
Factor <- list('CHL'=1, 'CAR'=1, 'EWT'=1000, 'LMA'=1000)
PlotObj <- list()
for (parm in names(Estimate_SpectralShifting)){
  Labs <- c(paste('Measured',parm,UnitsParms[[parm]]), 
            paste('Estimated',parm,UnitsParms[[parm]]))
  
  for (col0 in colnames(Estimate_SpectralShifting[[parm]])){
    print(c(parm,col0))
    valStep <- as.numeric(col0)
    if (valStep<10)
    {subChar <- paste('00',as.character(valStep),sep = '')}
    else if (valStep<100)
    {subChar <- paste('0',as.character(valStep),sep = '')}
    else if (valStep<1000)
    {subChar <- as.character(valStep)}
    fileName <- file.path(SpectralShifting_subDir,'FIGURES', paste(parm,'_',subChar,'.png',sep = ''))
    PlotObj <- scatter_inversion(target = Factor[[parm]]*Biochemistry[[parm]],
                                 estimate = Factor[[parm]]*Estimate_SpectralShifting[[parm]][[col0]],
                                 Colors = PlotCols[[parm]],
                                 Labs = Labs,
                                 fileName = fileName,
                                 categories = "RT_FS",
                                 PlotStats = T,
                                 MinMaxAxis = MinMax[[parm]],
                                 size = 1)
  }
}

