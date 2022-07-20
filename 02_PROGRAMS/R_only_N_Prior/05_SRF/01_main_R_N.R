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
source('../../Libraries/Lib_Gaussian_filter.R')
source('../../Libraries/Lib_Analysis_Inversion.R')
source('../../Libraries/Lib_Plots.R')
################################################################################
# input output directories
################################################################################
PathData <- '../../../01_DATA'
PathResults <- '../../../03_RESULTS/R_only_N_Prior'
PathSRF <- file.path(PathResults,'05_SRF')
dir.create(path = PathSRF,showWarnings = F,recursive = T)

################################################################################
# repository where data are stored
################################################################################
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
load(PathLOPdb)
lambda <- Reflectance$V1
Reflectance <- Reflectance[,-1]
Transmittance<- Transmittance[,-1]
nbSamples <- ncol(Reflectance)

################################################################################
# Sensor simulation
################################################################################

################################################################################
# Perform PROSPECT inversion using simulated spectral configuration
################################################################################
List_FWHM <- list(1,2,3,4,5,10,15,20)
# save results
DirSave_Root <- file.path(PathSRF,'SRF_fwhm_')

#component's wavelenghts
wavelenght<-list()
wavelenght$CHL <- c(714,849)
wavelenght$CAR <- c(809,731,549,523)
wavelenght$EWT <- c(1393)
wavelenght$LMA <- c(1885,1729)


Parms2Estimate<-list()
Parms2Estimate$CHL<-Parms2Estimate$CAR <- c("CHL", "CAR","ANT","EWT","LMA")
Parms2Estimate$EWT<-Parms2Estimate$LMA <- c("LMA","EWT")
# define parameters to estimate during inversion

N_prior <- Get_Nprior(SpecPROSPECT = SpecPROSPECT, lambda = lambda, Refl = Reflectance, Tran = NULL)

# define parameters to estimate during inversion
InitValues <- list()
InitValues<-list()
for (i in c(1:nbSamples)) {
  InitValues[[i]] <- data.frame("CHL"=40, "CAR"=10, "ANT"=0.1, "BROWN"=0, "EWT"=0.01, "LMA"=0.01, "N"=N_prior[i,])
}

# checking 
List_FWHM_to_do <- list()
for (fwhm in List_FWHM){
  SRF_Dir <- paste(DirSave_Root, as.character(fwhm), 'nm',sep = '')
  if (!file.exists(file.path(SRF_Dir,"Results_Invert_PROSPECT.RData"))){
    List_FWHM_to_do<- append(List_FWHM_to_do,fwhm)
  }
}

for (fwhm in List_FWHM_to_do) {
  for (parm in c("CHL","CAR","EWT","LMA")) {
    Invert_est<-list()
    FilePath <- file.path(PathSRF,paste(parm,'_SRF', sep = ''))
    dir.create(path = FilePath,showWarnings = F,recursive = T)
      minWL <- min(lambda)
      maxWL <- max(lambda)
      # define spectral sampling based on minWL, minWL_OUT, FWHM
      wl <- wavelenght                          #wavelength range
      SpectralProps <- data.frame('wl'= unlist(wl),'fwhm'=fwhm)#wl group by fwhm in data.frame
      SensorName <-paste(parm,'_Filter_',fwhm,sep='')
      SRF <- prosail::GetRadiometry(SpectralProps = SpectralProps,
                                    SensorName = SensorName,
                                    Path_SensorResponse = SRF_Dir,
                                    SaveSRF = FALSE)
      # Inversion
      SpecPROSPECT_Sensor <- data.frame(prosail::applySensorCharacteristics(wvl = SpecPROSPECT$lambda,
                                                                            InRefl = SpecPROSPECT,
                                                                            SRF = SRF))
      RT_filter <- gaussian_filter(SRF = SRF, 
                                   lambda = lambda, 
                                   Refl=Reflectance, 
                                   Tran=Transmittance)
      message(paste('FWHM = ',fwhm, parm))
      pb <- progress::progress_bar$new(
        format = "PROSPECT inversion [:bar] :percent in :elapsedfull",
        total = nbSamples, clear = FALSE, width= 100)
      for (i in 1:nbSamples) {                                         #filtered reflectance-transmittance list
        pb$tick()
        Rtmp <- RT_filter$Refl[,i]
        Ttmp <- NULL
        Invert <- Invert_PROSPECT(SpecPROSPECT =  SpecPROSPECT_Sensor,
                                      Refl = Rtmp,
                                      Tran = Ttmp,
                                      Parms2Estimate = Parms2Estimate[[parm]],
                                      PROSPECT_version = 'D',
                                      InitValues = InitValues[[i]],)
        Invert_est<-rbind(Invert_est[[parm]],Invert)
      }
      save(Invert_est, file = file.path(PathSRF,paste(parm,"_SRF", sep = ""),paste("fwhm_",fwhm,".RData", sep = "")))
    }
}

# 
# i <- 0
# for (fwhm in List_FWHM_to_do){
#   i <- i + 1
#   SRF_Dir <- paste(DirSave_Root, as.character(fwhm), 'nm',sep = '')
#   dir.create(path = SRF_Dir ,showWarnings = F,recursive = T)
#   path_save_data<-file.path(SRF_Dir,"Results_Invert_PROSPECT.RData")
#   PROSPECT_inversion <- Invert_est[[i]]
#   save(PROSPECT_inversion, file = path_save_data)
#   
#   path_save_plot<-file.path(SRF_Dir,"Results_Invert_PROSPECT.png")
#   
#   plots<-scatter_inversion(target = Biochemistry$CHLa+Biochemistry$CHLb,
#                            estimate = PROSPECT_inversion$CHL,
#                            Labs = c("Mesured CHL (µg/cm²)","Estimated CHL (µg/cm²)"),
#                            categories = T,
#                            Colors = "#66CC00",
#                            PlotStats = T,
#                            fileName = path_save_plot)
#   
# }





# #############################################################################################
# ###########################################################################################
# ################################################################################
# # results data paths
# ################################################################################
# pathdata1 <- "../../03_RESULTS/R_only_N_Prior/05_SRF/SRF_fwhm_"
# pathdata3 <- "nm/Results_Invert_PROSPECT.RData"
# 
# nrmse_df = data.frame("fwhm"=NA, "nrmse" = NA)
# 
# for (i in c(1,2,3,4,5,10,15,20)) {
#   pathdata2 <- as.character(i)
#   pathdata01 <- paste(pathdata1,pathdata2,sep = "")
#   pathdata02 <- paste(pathdata01,pathdata3, sep = "")
#   load(pathdata02)
#   
#   nrmse_df = rbind(nrmse_df,data.frame("fwhm"= i, 
#                                        "nrmse" = unlist(get_performances_inversion(target = (Biochemistry$CHLa+Biochemistry$CHLb),
#                                                                                    estimate = PROSPECT_inversion$CHL,
#                                                                                    categories = T)[3])))
#   
# }
# nrmse_df <- nrmse_df[-c(1),]
# 
# my.formula <- y~ poly(x, 2, raw = TRUE)
# 
# plot <- ggplot(nrmse_df, aes(x = fwhm, y = nrmse)) +
#   geom_point(aes(x = fwhm, y = nrmse), colour = "#66CC00", size = 2) +
#   labs(x = 'CHL_fwhm (nm)',y = 'NRMSE') +
#   theme(legend.position = "bottom",        axis.text = element_text(size = 12),
#         axis.title.x = element_text(size = 14, face = "bold"),
#         axis.title.y = element_text(size = 14, face = "bold"))  + 
#   theme_bw() +
#   geom_smooth(aes(x = fwhm, y = nrmse, color = "#66CC00"),method = lm, formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +
#   ggpmisc::stat_poly_eq(formula = my.formula, 
#                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                         parse = TRUE)+
#   ylim(0,70)
# filename = file.path(PathResults,'05_SRF/NRMSE_FWHM_standard.png')
# ggsave(filename,plot = plot, device = "png", path = NULL,
#        scale = 1, width = 20, height = 13, units = "cm",
#        dpi = 600)
# 
# plot
