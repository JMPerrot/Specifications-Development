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
source('../Libraries/Lib_Gaussian_filter.R')
source('../Libraries/Lib_Analysis_Inversion.R')
source('../Libraries/Lib_Plots.R')
################################################################################
# input output directories
################################################################################
PathData <- '../../01_DATA'
PathResults <- '../../03_RESULTS'
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
List_FWHM <- list(1,2,3,4,5,10,15,20)
# save results
DirSave_Root <- file.path(PathSRF,'SRF_fwhm_')

wavelenght_tested <- c(714,849)
List_FWHM_to_do <- list()
for (fwhm in List_FWHM){
  SRF_Dir <- paste(DirSave_Root, as.character(fwhm), 'nm',sep = '')
  if (!file.exists(file.path(SRF_Dir,"Results_Invert_PROSPECT.RData"))){
    List_FWHM_to_do<- append(List_FWHM_to_do,fwhm)
  }
}

PROSPECT_inv <- lapply(X = List_FWHM_to_do,FUN =  Invert_PROSPECT_GaussianFilter_fwhm,
                       Sampling_in = lambda, Refl = Refl, Tran = Tran,
                       wavelenght_tested= wavelenght_tested)

i <- 0
for (fwhm in List_FWHM_to_do){
  i <- i + 1
  SRF_Dir <- paste(DirSave_Root, as.character(fwhm), 'nm',sep = '')
  dir.create(path = SRF_Dir ,showWarnings = F,recursive = T)
  path_save_data<-file.path(SRF_Dir,"Results_Invert_PROSPECT.RData")
  PROSPECT_inversion <- PROSPECT_inv[[i]]
  save(PROSPECT_inversion, file = path_save_data)
  
  path_save_plot<-file.path(SRF_Dir,"Results_Invert_PROSPECT.png")
  
  plots<-scatter_inversion(target = Biochemistry$CHLa+Biochemistry$CHLb,
                          estimate = PROSPECT_inversion$CHL,
                          Labs = c("Mesured CHL (µg/cm²)","Estimated CHL (µg/cm²)"),
                          categories = T,
                          Colors = "#66CC00",
                          PlotStats = T,
                          fileName = path_save_plot)

}





#############################################################################################
###########################################################################################
################################################################################
# results data paths
################################################################################
pathdata1 <- "../../03_RESULTS/06_SRF/SRF_fwhm_"
pathdata3 <- "nm/Results_Invert_PROSPECT.RData"

nrmse_df = data.frame("fwhm"=NA, "nrmse" = NA)

for (i in c(1,2,3,4,5,10,15,20)) {
  pathdata2 <- as.character(i)
  pathdata01 <- paste(pathdata1,pathdata2,sep = "")
  pathdata02 <- paste(pathdata01,pathdata3, sep = "")
  load(pathdata02)
  
  nrmse_df = rbind(nrmse_df,data.frame("fwhm"= i, 
                                       "nrmse" = unlist(get_performances_inversion(target = (Biochemistry$CHLa+Biochemistry$CHLb),
                                                                            estimate = PROSPECT_inversion$CHL,
                                                                            categories = T)[3])))
  
}
nrmse_df <- nrmse_df[-c(1),]

my.formula <- y~ poly(x, 2, raw = TRUE)

plot <- ggplot(nrmse_df, aes(x = fwhm, y = nrmse)) +
  geom_point(aes(x = fwhm, y = nrmse), colour = "#66CC00", size = 2) +
  labs(x = 'CHL_fwhm (nm)',y = 'NRMSE') +
  theme(legend.position = "bottom",        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))  + 
  theme_bw() +
  geom_smooth(aes(x = fwhm, y = nrmse, color = "#66CC00"),method = lm, formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +
  ggpmisc::stat_poly_eq(formula = my.formula, 
                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                        parse = TRUE)+
  ylim(0,70)
filename = file.path('../../03_RESULTS/06_SRF/NRMSE_FWHM_standard.png')
ggsave(filename,plot = plot, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plot
