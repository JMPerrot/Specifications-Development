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
library(doFuture)
library(progress)
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
SFS_Dir <- file.path(PathResults,'04_FeatureSelection')
dir.create(path = SFS_Dir,showWarnings = F,recursive = T)
SpectralSampling_Dir <- file.path(PathResults,'02_SpectralSampling')

################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
lambda <- Reflectance$V1
Reflectance <- Reflectance[,-1]
Transmittance <- Transmittance[,-1]

################################################################################
# Which spectral sampling to be used?
################################################################################
# CHL & CAR: get spectral sampling corresponding to min NRMSE value 
# EWT & LMA: get trade-off = 45
Parms2Estimate <- c('CHL','CAR','EWT','LMA')
# Parms2Estimate <- c('LMA')
Opt_Sampling <- list()
Stats <- list()
for (parm in Parms2Estimate){
  FileName <- file.path(SpectralSampling_Dir,paste(parm,'_Statistics.csv',sep = ''))
  Stats[[parm]] <- readr::read_delim(file = FileName,delim = '\t')
  if (parm =='CHL' | parm =='CAR'){
    Opt_Sampling[[parm]] <- Stats[[parm]]$Sampling[which(Stats[[parm]]$NRMSE==min(Stats[[parm]]$NRMSE))]
  } else if (parm =='EWT' | parm =='LMA'){
    Opt_Sampling[[parm]] <- 45
  }
}

################################################################################
# Define parameters for inversion
################################################################################
# define parameters to estimate during inversion
ParmsEstInv <- list()
ParmsEstInv$CHL <- c('CHL', 'CAR', 'EWT', 'LMA', 'N')
ParmsEstInv$CAR <- c('CHL', 'CAR', 'EWT', 'LMA', 'N')
ParmsEstInv$EWT <- c('EWT', 'LMA', 'N')
ParmsEstInv$LMA <- c('EWT', 'LMA', 'N')

# define parameters to estimate during inversion
InitValues <- list()
InitValues$CHL <- data.frame(CHL=45, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
InitValues$CAR <- data.frame(CHL=45, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
InitValues$EWT <- data.frame(CHL=45, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
InitValues$LMA <- data.frame(CHL=45, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)

################################################################################
# Perform SFS reduce spectral information
################################################################################
# number of CPU available
nbWorkers <- 8

# define spectral domain
SpectralDomain <- list()
SpectralDomain$CHL <- SpectralDomain$CAR <- list('minWL' = 400, 'maxWL' = 900)
SpectralDomain$EWT <- SpectralDomain$LMA <- list('minWL' = 1300, 'maxWL' = 2400)

Evol_NRMSE <- DiscardedWL <- list()
for (parm in Parms2Estimate){
  ## Perform SFS on initial set of spectral bands listWL
  listWL <- as.list(seq(from = SpectralDomain[[parm]]$minWL,
                        to = SpectralDomain[[parm]]$maxWL,
                        by = Opt_Sampling[[parm]]))
  
  # number of wavelengths
  TotalWL <- length(listWL)
  # initialize selected features to none
  SelectedWL <- unlist(listWL)
  # number of features to explore during SFS (decreases at each iteration)
  AllWL <- SelectedWL
  DiscardedWL[[parm]] <- c()
  # initialize list of results to keep
  ResWL <- NRMSE_Est <- list()
  Evol_NRMSE[[parm]] <- c()
  
  # perform SFS using multiprocessing
  registerDoFuture()
  plan(multisession, workers = nbWorkers)
  pb <- progress_bar$new(
    format = "PROSPECT inversion [:bar] :percent in :elapsedfull",
    total = TotalWL, clear = FALSE, width= 100)
  for (nbvars2select in 1:(TotalWL-1)){
    pb$tick()
    nbVars <- length(AllWL)
    NumVar_list <- as.list(seq(1,length(AllWL)))
    ResWL0 <- as.list(seq(1,length(AllWL)))
    backward_SFS <- function() {
      foreach(numvar = NumVar_list) %dopar% {
        # eliminate spectral band
        SelectedWLTmp <- SelectedWL[-numvar]
        # adjust spectral information
        DataFit <- FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                                   lambda = lambda,
                                   Refl = Reflectance,
                                   Tran = Transmittance,
                                   UserDomain = SelectedWLTmp,
                                   UL_Bounds = FALSE)
        # perform PROSPECT inversion
        Invert_est <- prospect::Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                                Refl = DataFit$Refl,
                                                Tran = DataFit$Tran,
                                                PROSPECT_version = 'D',
                                                Parms2Estimate = ParmsEstInv[[parm]],
                                                InitValues = InitValues[[parm]],progressBar = FALSE)
        # compute performances for target parameter
        Stats_inversion <- get_performances_inversion(target = Biochemistry[[parm]],
                                                      estimate = Invert_est[[parm]], 
                                                      categories= TRUE)
        # NRMSE is the performance metric to minimize
        return(Stats_inversion)
      }
    }
    subSFS <- backward_SFS()
    # get list of NRMSE obtained when 
    NRMSE_Est[[nbvars2select]] <- matrix(unlist(lapply(subSFS,'[[',3)),ncol = 1)
    rownames(NRMSE_Est[[nbvars2select]]) <- AllWL
    # select wavelength resulting in minimum NRMSE when discarded
    SelVar <- which(NRMSE_Est[[nbvars2select]]==min(NRMSE_Est[[nbvars2select]],na.rm = T))
    WhichVar <- as.numeric(rownames(NRMSE_Est[[nbvars2select]])[SelVar[1]])
    DiscardedWL[[parm]] <- c(DiscardedWL[[parm]],SelectedWL[SelVar])
    # delete selected component from AllWL
    AllWL <- AllWL[-which(AllWL==SelectedWL[SelVar])]
    # remove selected wavelength
    SelectedWL <- SelectedWL[-SelVar]
    Evol_NRMSE[[parm]] <- c(Evol_NRMSE[[parm]],min(NRMSE_Est[[nbvars2select]],na.rm = T))
  }
  plan(sequential)
}

# perform inversion for last spectral band
for (parm in Parms2Estimate){
  ## Perform SFS on initial set of spectral bands listWL
  listWL <- seq(from = SpectralDomain[[parm]]$minWL,
                to = SpectralDomain[[parm]]$maxWL,
                by = Opt_Sampling[[parm]])
  Last_WL <- which(is.na(match(listWL,DiscardedWL[[parm]])))
  DiscardedWL[[parm]] <- c(DiscardedWL[[parm]],listWL[Last_WL])
  AllFeat <- Stats[[parm]]$NRMSE[which(Stats[[parm]]$Sampling==Opt_Sampling[[parm]])]
  Evol_NRMSE[[parm]] <- c(AllFeat,Evol_NRMSE[[parm]])
}


for (parm in Parms2Estimate){
  df <- data.frame('Discarded_WL' = DiscardedWL[[parm]],
                   'NRMSE' = Evol_NRMSE[[parm]])
  filename <- file.path(SFS_Dir,paste(parm,'_FeatureSelection.csv',sep = ''))
  write_delim(x = df,
              file = filename,
              delim = '\t',
              col_names = T)
}

# for (parm in Parms2Estimate){
#   filename <- file.path(SFS_Dir,paste(parm,'_FeatureSelection.csv',sep = ''))
#   df <- read_delim(file = filename,
#                    delim = '\t',
#                    col_names = T)
#   DiscardedWL[[parm]] <- df$Discarded_WL
#   Evol_NRMSE[[parm]] <- df$NRMSE
# }

#   
#   
#   
#   
# 
#   nbVars <- length(listWL)
#   NumVar_list <- as.list(seq(1,length(AllVars)))
#   ResVar0 <- as.list(seq(1,length(AllVars)))
#   
#   
#   subfeatures_SFS <- function() {
#     foreach(numvar = NumVar_list) %dopar% {
#       SelectedVarsTmp <- c(SelectedVars,AllVars[[numvar]])
#       ResVar0 <- lapply(X = Spectralpop,FUN = spectral_variance_subset,
#                         vars = SelectedVarsTmp, stand='None',
#                         var_crown = F,var_residu = F)
#       
#       ResVar0 <- as.data.frame(do.call(rbind,(lapply(ResVar0,unlist)))[,1:2])
#       Corr_tot0 <- cor.test(ResVar0$var_tot,
#                             All_Populations_diversite$Shannon)$estimate
#       NRMSE_Est0 <- cor.test(ResVar0$var_sp,
#                                 All_Populations_diversite$Shannon)$estimate
#       return(list('ResVar0'=ResVar0,'NRMSE_Est0'=NRMSE_Est0,'Corr_tot0'=Corr_tot0))
#     }
#   }
#   subSFS <- subfeatures_SFS()
#   
#   ResVar[[nbvars2select]] <- lapply(subSFS,'[[',1)
#   Corr_species[[nbvars2select]] <- matrix(unlist(lapply(subSFS,'[[',2)),ncol = 1)
#   Corr_tot[[nbvars2select]]  <- matrix(unlist(lapply(subSFS,'[[',3)),ncol = 1)
#   rownames(Corr_tot[[nbvars2select]]) <- rownames(Corr_species[[nbvars2select]]) <- AllVars
#   # select component showing highest correlation between var_tot and Shannon from population
#   SelVar <- which(Corr_tot[[nbvars2select]]==max(Corr_tot[[nbvars2select]],na.rm = T))
#   WhichVar <- as.numeric(rownames(Corr_tot[[nbvars2select]])[SelVar])
#   # add selected component to selected vars
#   SelectedVars <- c(SelectedVars,WhichVar)
#   # delete selected component from AllVars
#   AllVars <- AllVars[-which(AllVars==WhichVar)]
#   EvolCorr <- c(EvolCorr,max(Corr_tot[[nbvars2select]],na.rm = T))
#   
#   
# 
#   
#   
#     
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# parms_est_df=data.frame("CAR"=NA, "deleted_wl" = NA , "obs" = NA)              #           
# 
# seqWL<-c(500:540)#seq(400, 1000, by =10)
# 
# 
# while(length(seqWL)>5){
#   nrmse_df = data.frame("deleted_wl"=NA, "NRMSE" = NA)
#   parms_est_df=data.frame("CAR"=NA, "deleted_wl" = NA , "obs" = NA) 
#   for (i in c(1:length(seqWL))){
#     seqWLR<-seqWL[-i] #sequence de longueure d'onde reduite
#     lambda<-c(unlist(Reflectance[,1]))
#     WLselect<-match(seqWLR,lambda)
#     
#     pb <- progress::progress_bar$new(
#       format = "PROSPECT inversion [:bar] :percent in :elapsedfull",
#       total = length(Reflectance), clear = FALSE, width= 100)
#     
#     for (j in c(2:length(Reflectance))) {
#       DataFit<-FitSpectralData(SpecPROSPECT = SpecPROSPECT,
#                                lambda = SpecPROSPECT$lambda,
#                                Refl = matrix(unlist(Reflectance[[j]])),
#                                Tran = matrix(unlist(Transmittance[[j]])),
#                                UserDomain = seqWLR,
#                                UL_Bounds = FALSE)
#       
#       Invert_est<-Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
#                                   Refl = DataFit$Refl,
#                                   Tran = DataFit$Tran,
#                                   InitValues = data.frame(CHL = 40, CAR = 10, ANT = 0.1, BROWN = 0.01, EWT = 0.01, LMA
#                                                           = 0.01, PROT = 0.001, CBC = 0.009, N = 1.5, alpha = 40),
#                                   Parms2Estimate = c('CHL','CAR'))
#       
#       est = data.frame("CAR"=Invert_est$CAR, 
#                        "deleted_wl" = seqWL[i],
#                        "obs" = j)    
#       
#       parms_est_df=rbind(parms_est_df, est)
#       pb$tick()
#     }
#     parms<- parms_est_df %>% filter(parms_est_df$deleted_wl == seqWL[i])
#     nrmse_wl=get_performances_inversion(target = Biochemistry$CAR,
#                                         estimate = parms$CAR,
#                                         categories= TRUE)[3]
#     nrmse_df=rbind(nrmse_df, data.frame("deleted_wl" = seqWL[i],"NRMSE" = nrmse_wl[[1]]))
#   }
#   
#   nrmse_df<-nrmse_df[-1,]
#   min_nrmse<-min(as.matrix(nrmse_df$NRMSE))
#   nline<-match(1,match(nrmse_df$NRMSE,min_nrmse))
#   seqWL<-seqWL[-nline]
# }
# 
# nrmse_df<-nrmse_df[-1,]
