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
library(ggpubr)
source('../../Libraries/Lib_Plots.R')
source('../../Libraries/Lib_Analysis_Inversion.R')

################################################################################
# input output directories
################################################################################
PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/R_only_N_Prior'
SFS_Dir <- file.path(PathResults,'04_FeatureSelection')
dir.create(path = SFS_Dir,showWarnings = F,recursive = T)
SpectralSampling_Dir <- file.path(PathResults,'02_SpecSampling')
SpectralShifting_Dir <- file.path(PathResults,'03_SpecShifting')
Reference_Dir <- file.path(PathResults,'01_Reference')



## load leaf optics dataset ----------------------------------------------------

dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)


## compute statistics ----------------------------------------------------------

Parms2Estimate <- c('CHL','CAR','LMA','EWT')#
Stats_inversion_Ref <- Stats_inversion_SS <- list()
for (parm in Parms2Estimate){
  # load reference#1 for inversion
  FileName <- file.path(Reference_Dir,paste(parm,'_REFERENCE#1.RData',sep = ''))
  load(FileName)
  Ref1 <- ResultsInversion
  # load reference#2 for inversion
  FileName <- file.path(Reference_Dir,paste(parm,'_REFERENCE#2R.RData',sep = ''))
  load(FileName)
  Ref2 <- ResultsInversion
  
  # compute performances
  # ref#1
  Stats_inversion_Ref[[parm]] <- Stats_inversion_SS[[parm]] <- list()
  Stats_inversion_Ref[[parm]][['REF1']] <- get_performances_inversion(target = Ref1$measured,
                                                                      estimate = Ref1$estimated, 
                                                                      categories= TRUE)
  Stats_inversion_Ref[[parm]][['REF2']] <- get_performances_inversion(target = Ref2$measured,
                                                                      estimate = Ref2$estimated, 
                                                                      categories= TRUE)
}

################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
lambda <- Reflectance$V1
Reflectance <- Reflectance[,-1]
Transmittance <- Transmittance[,-1]
nbSamples<-ncol(Reflectance)

################################################################################
# Which spectral sampling to be used?
################################################################################
# CHL & CAR: get spectral sampling corresponding to min NRMSE value 
# EWT & LMA: get trade-off = 45


Opt_Sampling <- list()
Stats <- list()
for (parm in Parms2Estimate){
  FileName <- file.path(SpectralSampling_Dir,paste(parm,'_Statistics.csv',sep = ''))
  Stats[[parm]] <- readr::read_delim(file = FileName,delim = '\t')
  if (parm =='CHL' | parm =='CAR'){
    Opt_Sampling[[parm]] <- Stats[[parm]]$Sampling[which(Stats[[parm]]$NRMSE==min(Stats[[parm]]$NRMSE))]
  } else if (parm =='LMA'){
    Opt_Sampling[[parm]] <- 39#Stats[[parm]]$Sampling[which(Stats[[parm]]$NRMSE==min(Stats[[parm]]$NRMSE))]
  }else if (parm =='EWT'){
    Opt_Sampling[[parm]] <- 31#Stats[[parm]]$Sampling[which(Stats[[parm]]$NRMSE==min(Stats[[parm]]$NRMSE))]
  }
}

################################################################################
# Define parameters for inversion
################################################################################
# define parameters to estimate during inversion
ParmsEstInv <- list()
ParmsEstInv$CHL <- ParmsEstInv$CAR <- c('CHL', 'CAR','ANT', 'EWT', 'LMA')
ParmsEstInv$EWT <- ParmsEstInv$LMA <- c('EWT', 'LMA')
 
# define parameters to estimate during inversion

N_prior <- Get_Nprior(SpecPROSPECT = SpecPROSPECT, lambda = lambda, Refl = Reflectance, Tran = NULL)
# define parameters to estimate during inversion
InitValues <- list()
InitValues<-list()
for (i in c(1:nbSamples)) {
  InitValues[[i]] <- data.frame("CHL"=40, "CAR"=10, "ANT"=0.1, "BROWN"=0, "EWT"=0.01, "LMA"=0.01, "N"=N_prior[i,])
}

################################################################################
# Perform SFS reduce spectral information
################################################################################
# number of CPU available
nbWorkers <- 2

# define spectral domain
Opt_Shifting <- list()
for (parm in Parms2Estimate){
  FileName <- file.path(SpectralShifting_Dir,paste(parm,'_Statistics.csv',sep = ''))
  Stats[[parm]] <- readr::read_delim(file = FileName,delim = '\t')
  Opt_Shifting[[parm]] <- Stats[[parm]]$Sampling[which(Stats[[parm]]$NRMSE==min(Stats[[parm]]$NRMSE))]
}
SpectralDomain <- list()
SpectralDomain$CHL <- list('minWL' = 400+Opt_Shifting$CHL, 'maxWL' = 900)
SpectralDomain$CAR <- list('minWL' = 400+Opt_Shifting$CAR, 'maxWL' = 900)
SpectralDomain$EWT <- list('minWL' = 1300+Opt_Shifting$EWT, 'maxWL' = 2400)
SpectralDomain$LMA <- list('minWL' = 1300+Opt_Shifting$LMA, 'maxWL' = 2400)

Evol_NRMSE <- DiscardedWL <- list()

PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
plotparm <- list()


for (parm in Parms2Estimate){
  filename <- file.path(SFS_Dir,paste(parm,'_FeatureSelection.csv',sep = ''))
  if (!file.exists(filename)){
    
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
    for (nbvars2select in c(1:(TotalWL-1))){
      pb$tick()
      nbVars <- length(AllWL)
      NumVar_list <- as.list(seq(1,length(AllWL)))
      ResWL0 <- as.list(seq(1,length(AllWL)))
      backward_SFS <-function() {
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
          Invert_est<-list()
          for (i in c(1:nbSamples)) {
            Invert <- prospect::Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                                Refl = DataFit$Refl[[i]],
                                                Tran = NULL,
                                                PROSPECT_version = 'D',
                                                Parms2Estimate = ParmsEstInv[[parm]],
                                                InitValues = InitValues[[i]])
            Invert_est<-rbind(Invert_est,Invert)
          }
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
    
    
    # perform inversion for last spectral band
    
    ## Perform SFS on initial set of spectral bands listWL
    listWL <- seq(from = SpectralDomain[[parm]]$minWL,
                  to = SpectralDomain[[parm]]$maxWL,
                  by = Opt_Sampling[[parm]])
    Last_WL <- which(is.na(match(listWL,DiscardedWL[[parm]])))
    DiscardedWL[[parm]] <- c(DiscardedWL[[parm]],listWL[Last_WL])
    AllFeat <- Stats[[parm]]$NRMSE[which(Stats[[parm]]$Sampling==Opt_Sampling[[parm]])]
    Evol_NRMSE[[parm]] <- c(AllFeat,Evol_NRMSE[[parm]])
    
    
    
    
    SpecSampling <- data.frame('Discarded_WL' = DiscardedWL[[parm]],
                     'NRMSE' = Evol_NRMSE[[parm]])
    
    
    write_delim(x = SpecSampling,
                file = filename,
                delim = '\t',
                col_names = T)
  }else{
    SpecSampling <- readr::read_delim(file = filename,delim = '\t')
  }
  # load inversion results for spectral samplings
  
  fileplot<-file.path(SFS_Dir,paste(parm,'_FeatureSelection.png',sep = ''))
  
  plotparm[[parm]]<-ggplot(SpecSampling, aes(x = seq(length(Discarded_WL),1, -1), y = NRMSE))+
    geom_line(aes(x = seq(length(Discarded_WL),1, -1), y = NRMSE), colour = PlotCols[[parm]], size = 1)+
    labs(x="number of wl selected",y="NRMSE (%)") +
    ylim(0,100)+
    theme_bw() +
    geom_abline(slope = 0, intercept = Stats_inversion_Ref[[parm]]$REF1$NRMSE, linetype='dashed',size=1,col='black') +
    geom_abline(slope = 0, intercept = Stats_inversion_Ref[[parm]]$REF2$NRMSE,linetype='dashed',size=1,col='green')
  
  filename_plot<- file.path(SFS_Dir,paste(parm,'_FeatureSelection.png',sep = ''))
  
  ggsave(filename_plot,plot = plotparm[[parm]], device = "png", path = NULL,
         scale = 1, width = 20, height = 13, units = "cm",
         dpi = 600)
}
plot_parm1 <- ggarrange(plotparm$CHL,plotparm$CAR,plotparm$EWT,plotparm$LMA,
                         plotlist = NULL,
                         ncol = 2,
                         nrow = 2,
                         labels = NULL,
                         label.x = 0,
                         label.y = 1,
                         hjust = -0.5,
                         vjust = 1.5,
                         font.label = list(size = 10, 
                                           color = "black", 
                                           face = "bold", 
                                           family = NULL),
                         align = "none",
                         widths = 100,
                         heights = 100,
                         legend = "none",
                         common.legend = TRUE
)
plot_parm<- ggpubr::annotate_figure(plot_parm1, 
                                     bottom = text_grob("NRMSE en du nombre de longueur d'onde", 
                                                        color = "black", 
                                                        face = "bold", 
                                                        size = 14))
## save figures ----------------------------------------------------------------
ggsave(file.path(SFS_Dir,
                 paste('NRMSE_ALL_','SFS.png',sep = '')), 
       plot_parm,
       device = "png")

