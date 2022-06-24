################################################################################
## this program allows to run the PROSPECT model by modifying only the FWHM 
## of the Gaussian filter. This allows to simulate sources with different 
## spectral widths in order to determine the maximum spectral width for 
## a prototype 
################################################################################
# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

################################################################################
## Libraries
################################################################################
library(prospect)
library(prosail)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(Metrics)
library(progress)
library(ggpubr)
source('../Libraries/Lib_Gaussian_filter.R')
source('../Libraries/Lib_Analysis_Inversion.R')
source('../Libraries/Lib_Plots.R')
################################################################################
# input output directories
################################################################################
PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/R&T'
PathSRF <- file.path(PathResults,'06_SRF')
dir.create(path = PathSRF,showWarnings = F,recursive = T)

################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
load(file.path(PathResults,'01_Reference/REF#1.RData'))

################################################################################
# Sensor simulation
################################################################################
lambda <- c(matrix(unlist(Reflectance[,1]),ncol = 1))
Refl <- matrix(unlist(Reflectance[,-1], use.names=FALSE),nrow = length(lambda))
Tran <- matrix(unlist(Transmittance[,-1], use.names=FALSE),nrow = length(lambda))

################################################################################
# Perform PROSPECT inversion using simulated spectral configuration
################################################################################
List_FWHM <- list(1,2,3,4,5,10,15,20,25)
# save results
DirSave_Root <- file.path(PathSRF,'SRF_fwhm_')
wavelenght_tested <- list()
wavelenght_tested$CHL <- c(714,849)
wavelenght_tested$CAR <- c(557,518)
wavelenght_tested$EWT <- c(1792,2332,2017,1882)
wavelenght_tested$LMA <- c(1882,1747,1477,2287)
List_FWHM_to_do <- list()
parm2estimate = c("CHL","CAR","LMA","EWT")
PROSPECT_inv<- list()

for (parm in parm2estimate) {
  for (fwhm in List_FWHM){
    SRF_Dir <- paste(DirSave_Root, as.character(fwhm), 'nm',sep = '')
    if (!file.exists(file.path(SRF_Dir,paste(parm,"_Results_Invert_PROSPECT.RData", sep = "")))){
      List_FWHM_to_do<- append(List_FWHM_to_do,fwhm)
    }else{
     load(file.path(SRF_Dir,paste(parm,"_Results_Invert_PROSPECT.RData", sep = "")))
    }
  }
  List_FWHM_to_do <- unlist(List_FWHM_to_do)
  # PROSPECT_inv[[parm]] <- lapply(X = List_FWHM_to_do,FUN =  Invert_PROSPECT_GaussianFilter_fwhm,
  #                      Sampling_in = lambda, Refl = Refl, Tran = Tran,
  #                      wavelenght_tested= wavelenght_tested[parm], parm = parm)

  i <- 0
  plots<-list()
  PlotCols <- list('CHL' = "#66CC00", 'CAR' = "orange", 'LMA' = "red", 'EWT' = "blue")
  for (fwhm in List_FWHM_to_do){
    PROSPECT_inv[[parm]]<-Invert_PROSPECT_GaussianFilter_fwhm(Sampling_in = lambda, Refl = Refl, Tran = Tran, fwhm = fwhm,
                                                             wavelenght_tested= wavelenght_tested[parm], parm = parm)
    i <- i + 1
    SRF_Dir <- paste(DirSave_Root, as.character(fwhm), 'nm',sep = '')
    dir.create(path = SRF_Dir ,showWarnings = F,recursive = T)
    path_save_data<-file.path(SRF_Dir,paste(parm,"_Results_Invert_PROSPECT.RData", sep = ""))
    PROSPECT_inversion <- PROSPECT_inv[[parm]]
    save(PROSPECT_inversion, file = path_save_data)
    
    path_save_plot<-file.path(SRF_Dir,paste(parm,"_Results_Invert_PROSPECT.png",spe = ""))
   # for(parm in c("CHL","CAR","EWT","LMA")){
      plots[[parm]]<-scatter_inversion(target = Biochemistry[[parm]],
                              estimate = PROSPECT_inversion[[parm]],
                              Labs = c(paste("Mesured", parm, "(µg/cm²)"),paste("Estimated", parm, "(µg/cm²)")),
                              categories = T,
                              Colors = PlotCols[[parm]],
                              PlotStats = T,
                              fileName = path_save_plot,)
    #}
  }
}




#############################################################################################
## compute statistics ----------------------------------------------------------
Reference_Dir <- file.path(PathResults,'01_Reference')
Parms2Estimate <- c('CHL','CAR','EWT','LMA')#
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
###########################################################################################
################################################################################
# results data paths
################################################################################
pathdata1 <- "../../../03_RESULTS/R&T/06_SRF/SRF_fwhm_"
pathdata3 <- "nm/"

plot<-list()
for(parm in parm2estimate){
  nrmse_df = data.frame("fwhm"=NA, "nrmse" = NA)
  for (i in List_FWHM){

    pathdata2 <- as.character(i)
    pathdata01 <- paste(pathdata1,pathdata2,sep = "")
    pathdata02 <- paste(pathdata01,pathdata3, sep = "")
    pathdata_final <- paste(pathdata02, parm, "_Results_Invert_PROSPECT.RData", sep = "")
    load(pathdata_final)
    nrmse_df = rbind(nrmse_df,data.frame("fwhm"= i,
                                         "nrmse" = unlist(get_performances_inversion(target = (Biochemistry[[parm]]),
                                                                              estimate = PROSPECT_inversion[[parm]],
                                                                              categories = T)[3])))
  }

  nrmse_df <- nrmse_df[-c(1),]
  
  my.formula <- y~ poly(x, 2, raw = TRUE)
  
  plot[[parm]] <- ggplot(nrmse_df, aes(x = fwhm, y = nrmse)) +
    geom_point(aes(x = fwhm, y = nrmse), colour = PlotCols[[parm]], size = 2) +
    labs(x = paste(parm, "fwlm (nm)") ,y = 'NRMSE (%)') +
    theme(legend.position = "bottom",        axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"))  +
    theme_bw() +
    geom_abline(slope = 0, intercept = Stats_inversion_Ref[[parm]]$REF1$NRMSE,
                linetype='dashed',size=1,col='black') +
    geom_abline(slope = 0, intercept = Stats_inversion_Ref[[parm]]$REF2$NRMSE,
                linetype='dashed',size=1,col='green')+
    geom_smooth(aes(x = fwhm, y = nrmse, color = "#66CC00"),method = lm, formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +
    ggpmisc::stat_poly_eq(formula = my.formula,
                          aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                          parse = TRUE)+
    ylim(0,90)
  
  filename = file.path('../../../03_RESULTS/R&T/06_SRF/', paste(parm, '_NRMSE_FWHM_standard.png', sep = ""))
  ggsave(filename,plot = plot[[parm]], device = "png", path = NULL,
         scale = 1, width = 20, height = 13, units = "cm",
         dpi = 600)
  
  plot
}
plot_nrmse1 <- ggarrange(plot$CHL,plot$CAR,plot$EWT,plot$LMA,
                         plotlist = NULL,
                         ncol = 4,
                         nrow = 1,
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
                                     bottom = text_grob("NRMSE en fonction de FWHM", 
                                                        color = "black", 
                                                        face = "bold", 
                                                        size = 14))
## save figures ----------------------------------------------------------------
ggsave(file.path('../../../03_RESULTS/R&T/06_SRF/All_NRMSE_FWHM_standard.png'), 
       plot_nrmse,
       device = "png")

