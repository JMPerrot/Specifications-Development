################################################################################
## This code summarizes the evolution of the performances of PROSPECT inversion
################################################################################
# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


## Libraries required ----------------------------------------------------------

library(tidyverse)
library(prospect)
library(data.table)
library(doFuture)
library(ggplot2)
library(ggpubr)

source('../Libraries/Lib_Plots.R')
source('../Libraries/Lib_Analysis_Inversion.R')


## Input output directories ----------------------------------------------------

PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/R_only'
Reference_Dir <- file.path(PathResults,'01_Reference')
SpectralSampling_Dir <- file.path(PathResults,'02_SpectralSampling')


## load leaf optics dataset ----------------------------------------------------

dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)


## compute statistics ----------------------------------------------------------

Parms2Estimate <- c('CHL','CAR','EWT','LMA')
Stats_inversion_Ref <- Stats_inversion_SS <- list()
for (parm in Parms2Estimate){
  # load reference#1 for inversion
  FileName <- file.path(Reference_Dir,paste(parm,'_REFERENCE#1.RData',sep = ''))
  load(FileName)
  Ref1 <- ResultsInversion
  # load reference#2 for inversion
  FileName <- file.path(Reference_Dir,paste(parm,'_REFERENCE#2.RData',sep = ''))
  load(FileName)
  Ref2 <- ResultsInversion
  # load inversion results for spectral samplings
  FileName <- file.path(SpectralSampling_Dir,paste(parm,'_SpecSampling.csv',sep = ''))
  SpecSampling <- readr::read_delim(file = FileName,delim = '\t')
  
  # compute performances
  # ref#1
  Stats_inversion_Ref[[parm]] <- Stats_inversion_SS[[parm]] <- list()
  Stats_inversion_Ref[[parm]][['REF1']] <- get_performances_inversion(target = Ref1$measured,
                                                              estimate = Ref1$estimated, 
                                                              categories= TRUE)
  Stats_inversion_Ref[[parm]][['REF2']] <- get_performances_inversion(target = Ref2$measured,
                                                              estimate = Ref2$estimated, 
                                                              categories= TRUE)
  for (ss in colnames(SpecSampling)){
    Stats_inversion_SS[[parm]][[ss]] <- get_performances_inversion(target = Biochemistry[[parm]],
                                                                   estimate = SpecSampling[[ss]],
                                                                   categories= TRUE)
  }
  Stats_inversion_SS[[parm]] <- do.call(rbind,Stats_inversion_SS[[parm]])
  Stats_inversion_SS[[parm]]$Sampling <- as.numeric(rownames(Stats_inversion_SS[[parm]]))
}

## save Statistics -------------------------------------------------------------
for (parm in names(Stats_inversion_SS)){
  FileName <- file.path(SpectralSampling_Dir,paste(parm,'_Statistics.csv',sep = ''))
  write_delim(x = Stats_inversion_SS[[parm]],
              file = FileName,
              delim = '\t',
              col_names = T)
}



## produce figure --------------------------------------------------------------


PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")

plot0 <- list()
for (parm in Parms2Estimate){
  plot0[[parm]] <- ggplot(Stats_inversion_SS[[parm]], aes(x = Sampling, y = NRMSE)) +
    geom_line(aes(x = Sampling, y = NRMSE), colour = PlotCols[[parm]], size = 1) +
    labs(x="Spectral sampling (nm)",y="NRMSE (%)") +
    ylim(0, max(Stats_inversion_SS[[parm]]$NRMSE)) +
    theme_bw() +
    geom_abline(slope = 0, intercept = Stats_inversion_Ref[[parm]]$REF1$NRMSE,
                linetype='dashed',size=1,col='black') +
    geom_abline(slope = 0, intercept = Stats_inversion_Ref[[parm]]$REF2$NRMSE,
                linetype='dashed',size=1,col='green')
  
  filename = file.path(SpectralSampling_Dir,paste('NRMSE_',parm,'_SpectralSampling.png',sep = ''))
  ggsave(filename,plot = plot0[[parm]], device = "png", path = NULL,
         scale = 1, width = 20, height = 13, units = "cm",
         dpi = 600)
}
plot_nrmse1 <- ggarrange(plot0$CHL,plot0$CAR,plot0$EWT,plot0$LMA,
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
plot_nrmse<- ggpubr::annotate_figure(plot_nrmse1, 
                             bottom = text_grob("NRMSE en fonction du pas", 
                                                color = "black", 
                                                face = "bold", 
                                                size = 14))

## save figures ----------------------------------------------------------------
ggsave(file.path(SpectralSampling_Dir,
                 paste('NRMSE_ALL','_SpectralSampling.png',sep = '')), 
       plot_nrmse,
       device = "png")
