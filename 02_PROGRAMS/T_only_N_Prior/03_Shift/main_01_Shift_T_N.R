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
library(prosail)
source('../../Libraries/Lib_Plots.R')
source('../../Libraries/Lib_Analysis_Inversion.R')


# input output directories -----------------------------------------------------

PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/T_only_N_Prior'
SpecShifting_Dir <- file.path(PathResults,'03_SpecShifting')
SpecSampling_Dir <- file.path(PathResults,'02_SpecSampling')
dir.create(path = SpecShifting_Dir, showWarnings = F,recursive = T)







# load leaf optics dataset -----------------------------------------------------
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
lambda <- Reflectance$V1
Reflectance <- Reflectance[,-1]
Transmittance<- Transmittance[,-1]
nbSamples <- ncol(Reflectance)

# Define parameters for inversion ----------------------------------------------

Parms2Estimate <- c('LMA','EWT','CHL','CAR')
Parms2Estimate_ind <- list('CHL'='CHL','EWT'= 'EWT','LMA'='LMA','CAR'='CAR')

Stats<-list()
Opt_Sampling<-list()
for (parm in Parms2Estimate){
  FileName <- file.path(SpecSampling_Dir,paste(parm,'_Statistics.csv',sep = ''))
  Stats[[parm]] <- readr::read_delim(file = FileName,delim = '\t')
  Opt_Sampling[[parm]] <- Stats[[parm]]$Sampling[which(Stats[[parm]]$NRMSE==min(Stats[[parm]]$NRMSE))]
}
# define spectral sampling
SpecSampling <- list()
SpecSampling$CHL <- Opt_Sampling$CHL
SpecSampling$CAR <- Opt_Sampling$CAR
SpecSampling$LMA <- 39#Opt_Sampling$LMA
SpecSampling$EWT <- 31#Opt_Sampling$EWT
# define spectral domain
minlambda <- list()
minlambda$CHL <- c(0:SpecSampling$CHL)#max(unlist(SpecSampling)))
minlambda$EWT <- c(0:SpecSampling$EWT)#max(unlist(SpecSampling)))
minlambda$CAR <- c(0:SpecSampling$CAR)#max(unlist(SpecSampling)))
minlambda$LMA <- c(0:SpecSampling$LMA)#max(unlist(SpecSampling)))


# define parameters to estimate during inversion
ParmsEstInv <- list()
ParmsEstInv$CHL<-ParmsEstInv$CAR <- "ALL"
ParmsEstInv$EWT<-ParmsEstInv$LMA <- c('EWT', 'LMA')

# define parameters to estimate during inversion

N_prior <- Get_Nprior(SpecPROSPECT = SpecPROSPECT, lambda = lambda, Refl = NULL, Tran = Transmittance)


# define parameters to estimate during inversion
InitValues <- list()
InitValues<-list()
for (i in c(1:nbSamples)) {
  InitValues[[i]] <- data.frame("CHL"=40, "CAR"=10, "ANT"=0.1, "BROWN"=0, "EWT"=0.01, "LMA"=0.01, "N"=N_prior[i,])
}

# Perform inversion ------------------------------------------------------------

# for each parameter
nbWorkers <- 2
registerDoFuture()
plan(multisession, workers = nbWorkers)
Estimate_SpecShifting <- list()
SpecShifting_subDir <- file.path(SpecShifting_Dir,'_individualSteps')
dir.create(path = SpecShifting_subDir,
           showWarnings = F,recursive = T)

for (parm in Parms2Estimate){
  print(parm)
  # apply multiprocessing for each spectral sampling
  Invert_SpecShifting <- function() {
    foreach(subshifting = c(0:max(unlist(SpecSampling[[parm]])))) %dopar% {
      SpectralDomain <- list()
      SpectralDomain$CHL <- SpectralDomain$CAR <- list('minWL' = 400 +subshifting, 'maxWL' = 900) 
      SpectralDomain$EWT <- SpectralDomain$LMA <- list('minWL' = 1300+subshifting, 'maxWL' = 2400)
      if (subshifting<10){
        subChar <- paste('00',as.character(subshifting),sep = '')
      }else if (subshifting<100){
        subChar <- paste('0',as.character(subshifting),sep = '')
      }else if (subshifting<1000){
          subChar <- as.character(subshifting)
          }
      FileName <- file.path(SpecShifting_subDir,
                            paste(parm,'_SpecShifting_',subChar,'.RData',sep = ''))
      if (!file.exists(FileName)){
        lambda_tmp <- seq(from = SpectralDomain[[parm]]$minWL,
                          to = SpectralDomain[[parm]]$maxWL,
                          by = SpecSampling[[parm]])
        
        DataFit <- FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                                   lambda = unlist(lambda),
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
  PROSPECT_EST <- Invert_SpecShifting()
  for (subparm in Parms2Estimate_ind[[parm]]){
    Estimate_SpecShifting[[subparm]] <- do.call(cbind,lapply(PROSPECT_EST,'[[',subparm))
    Estimate_SpecShifting[[subparm]] <- data.frame(Estimate_SpecShifting[[subparm]])
    colnames(Estimate_SpecShifting[[subparm]]) <- minlambda[[parm]]
  }
}
plan(sequential)

# save all spectral sampling results in a unique file per parameter
for (parm in names(Estimate_SpecShifting)){
  FileName <- file.path(SpecShifting_Dir,paste(parm,'_SpecShifting.csv',sep = ''))
  write_delim(x = Estimate_SpecShifting[[parm]],
              file = FileName,
              delim = '\t',
              col_names = T)
}


## Produce figures ------------------------------------------------------------


# plot results
PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
MinMax <- list('CHL' = c(0,110), 'CAR'  = c(0,130), 'LMA' = c(0,165), 'EWT' = c(0,160))
UnitsParms <- list('CHL'='(µg/cm²)', 'CAR'='(µg/cm²)', 'EWT'='(mg/cm²)', 'LMA'='(mg/cm²)')
Factor <- list('CHL'= 1, 'CAR' = 1, 'EWT'= 1000, 'LMA'= 1000)
PlotObj <- list()
for (parm in names(Estimate_SpecShifting)){
  Labs <- c(paste('Measured',parm,UnitsParms[[parm]]), 
            paste('Estimated',parm,UnitsParms[[parm]]))
  
  for (col0 in colnames(Estimate_SpecShifting[[parm]])){
    print(c(parm,col0))
    valStep <- as.numeric(col0)
    if (valStep<10)
    {subChar <- paste('00',as.character(valStep),sep = '')}
    else if (valStep<100)
    {subChar <- paste('0',as.character(valStep),sep = '')}
    else if (valStep<1000)
    {subChar <- as.character(valStep)}
    fileName <- file.path(SpecShifting_subDir,'FIGURES', paste(parm,'_',subChar,'.png',sep = ''))
    PlotObj <- scatter_inversion(target = Factor[[parm]]*Biochemistry[[parm]],
                                 estimate = Factor[[parm]]*Estimate_SpecShifting[[parm]][[col0]],
                                 Colors = PlotCols[[parm]],
                                 Labs = Labs,
                                 fileName = fileName,
                                 categories = "T_N",
                                 PlotStats = T,
                                 MinMaxAxis = MinMax[[parm]],
                                 size = 1)
  }
}

