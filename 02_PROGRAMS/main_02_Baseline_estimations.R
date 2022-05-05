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
source('../Libraries/Lib_Analysis_Inversion.R')
source('../Libraries/Lib_Plots.R')

# repository where data are stored
gitlab_Rep <- 'https://gitlab.com/jbferet/myshareddata/raw/master/LOP/'

# download ANGERS data
dbName <- 'ANGERS'

# files available
fileName <- list('DataBioch.txt','ReflectanceData.txt','TransmittanceData.txt')
Biochemistry <- Refl <- Tran <- list()
Biochemistry <- fread(paste(gitlab_Rep,dbName,'/',fileName[[1]],sep=''))
Refl<- fread(paste(gitlab_Rep,dbName,'/',fileName[[2]],sep=''))
Tran <- fread(paste(gitlab_Rep,dbName,'/',fileName[[3]],sep=''))

# Get the wavelengths corresponding to the reflectance and transmittance measurements  
lambda <- unlist(Refl[,1], use.names=FALSE)
Refl <- matrix(unlist(Refl[,-1], use.names=FALSE),nrow = length(lambda))
Tran <- matrix(unlist(Tran[,-1], use.names=FALSE),nrow = length(lambda))

# Get the number of samples
nbSamples <- ncol(Refl)

# Estimate all parameters for PROSPECT-D
Parms2Estimate  = 'ALL'
CHL_ALL <- CAR_ALL <- ANT_ALL <- EWT_ALL <- LMA_ALL <- N_ALL <- c()
InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
# Adjust spectral domain for SpecPROSPECT to fit leaf optical properties 
SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl, 
                           Tran =Tran,UserDomain = c(lambda[1],lambda[length(lambda)]),UL_Bounds = TRUE)
SubSpecPROSPECT = SubData$SpecPROSPECT
Sublambda       = SubData$lambda
SubRefl         = SubData$Refl
SubTran         = SubData$Tran

print('PROSPECT inversion using full spectral range')
for (i in 1:nbSamples[[1]]){
  print(i)
  res <- Invert_PROSPECT(SubSpecPROSPECT, Refl = SubRefl[,i], Tran = SubTran[,i], 
                         PROSPECT_version = 'D',Parms2Estimate = Parms2Estimate, 
                         InitValues = InitValues)
  CHL_ALL[i] <- res$CHL
  CAR_ALL[i] <- res$CAR
  ANT_ALL[i] <- res$ANT
  EWT_ALL[i] <- res$EWT
  LMA_ALL[i] <- res$LMA
  N_ALL[i] <- res$N
}

NRMSE<- get_performances_inversion(target = Biochemistry$CHLa+Biochemistry$CHLb,
                                   estimate = CHL_ALL,
                                   categories = T)[[3]]

X<-data.frame(Biochemistry$CHLa+Biochemistry$CHLb,CHL_ALL)
save(res,file = "../../03_RESULTS/REFERENCE/REF#1.RData")

CHL_1 = scatter_inversion(target = Biochemistry$CHLa+Biochemistry$CHLb,
                          estimate = CHL_ALL,
                          Colors = "#66CC00",
                          Labs = c("Mesured CHL (µg/cm²)","Estimated CHL (µg/cm²)"),
                          fileName = "../../03_RESULTS/REFERENCE/FIG/CHL_REF#1.png",
                          categories = "RT_FS",
                          PlotStats = T)
CAR_1 = scatter_inversion(target = Biochemistry$CAR,
                          estimate = CAR_ALL,
                          Colors = "orange",
                          Labs = c("Mesured CAR (µg/cm²)","Estimated CAR (µg/cm²)"),
                          fileName = "../../03_RESULTS/REFERENCE/FIG/CAR_REF#1.png",
                          categories = "RT_FS",
                          PlotStats = T)
LMA_1 = scatter_inversion(target = Biochemistry$LMA,
                          estimate = LMA_ALL,
                          Colors = "red",
                          Labs = c("Mesured LMA (g/cm²)","Estimated LMA (g/cm²)"),
                          fileName = "../../03_RESULTS/REFERENCE/FIG/LMA_REF#1.png",
                          categories = "RT_FS",
                          PlotStats = T)
EWT_1 = scatter_inversion(target = Biochemistry$EWT,
                          estimate = EWT_ALL,
                          Colors = "blue",
                          Labs = c("Mesured EWT (g/cm²)","Estimated EWT (g/cm²)"),
                          fileName = "../../03_RESULTS/REFERENCE/FIG/EWT_REF#1.png",
                          categories = "RT_FS",
                          PlotStats = T)

plotREFa1<-ggarrange(CHL_1,CAR_1,EWT_1,LMA_1,
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
plotREF1<- annotate_figure(plotREFa1, 
                           bottom = text_grob("Données de référence REF#1", 
                                              color = "black", 
                                              face = "bold", 
                                              size = 14))
ggsave("../../03_RESULTS/REFERENCE/REF#1.png", plotREF1,device = "png")
################################################################################

# Estimate all parameters for PROSPECT-D
Parms2Estimate  = c('CHL','CAR','ANT','EWT','LMA')
# Parms2Estimate  = c('LMA')
InitValues <- data.frame(CHL=40, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
print('PROSPECT inversion using optimal setting')
ParmEst <- Invert_PROSPECT_OPT(SpecPROSPECT, lambda=lambda, Refl = Refl,
                               Tran = Tran, PROSPECT_version = 'D',
                               Parms2Estimate = Parms2Estimate, InitValues = InitValues)
CHL_OPT <- ParmEst$CHL
CAR_OPT <- ParmEst$CAR
ANT_OPT <- ParmEst$ANT
EWT_OPT <- ParmEst$EWT
LMA_OPT <- ParmEst$LMA


save(ParmEst,file = "../../03_RESULTS/REFERENCE/REF#2.RData")

CHL_2 = scatter_inversion(target = Biochemistry$CHLa+Biochemistry$CHLb,
                          estimate = CHL_OPT,
                          Colors = "#66CC00",
                          fileName = "../../03_RESULTS/REFERENCE/FIG/CHL_REF#2.png",
                          Labs = list("Mesured CHL (µg/cm²)","Estimated CHL (µg/cm²)"),
                          PlotStats = TRUE,
                          categories = "RT_OPT")
CAR_2 = scatter_inversion(target = Biochemistry$CAR,
                          estimate = CAR_OPT,
                          Colors = "orange",
                          Labs = c("Mesured CAR (µg/cm²)","Estimated CAR (µg/cm²)"),
                          fileName = "../../03_RESULTS/REFERENCE/FIG/CAR_REF#2.png",
                          categories = "RT_OPT",
                          PlotStats = T)
LMA_2 = scatter_inversion(target = Biochemistry$LMA,
                          estimate = LMA_OPT,
                          Colors = "red",
                          Labs = c("Mesured LMA (g/cm²)","Estimated LMA (g/cm²)"),
                          fileName = "../../03_RESULTS/REFERENCE/FIG/LMA_REF#2.png",
                          categories = "RT_OPT",
                          PlotStats = T)
EWT_2 = scatter_inversion(target = Biochemistry$EWT,
                          estimate = EWT_OPT,
                          Colors = "blue",
                          Labs = c("Mesured EWT (g/cm²)","Estimated EWT (g/cm²)"),
                          fileName = "../../03_RESULTS/REFERENCE/FIG/EWT_REF#2.png",
                          categories = "RT_OPT",
                          PlotStats = T)

plotREFa2 <- ggarrange(CHL_2,CAR_2,EWT_2,LMA_2,
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
plotREF2<- annotate_figure(plotREFa2, 
                           bottom = text_grob("Données de référence REF#2", 
                                              color = "black", 
                                              face = "bold", 
                                              size = 14))

ggsave("../../03_RESULTS/REFERENCE/REF#2.png", plotREF2,device = "png")
