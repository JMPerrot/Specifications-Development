################################################################################
## this code calculates the RMSE and displays it for the different shift 
## (offset) values
################################################################################

# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
################################################################################
# input function and data 
################################################################################
library(ggplot2)
library(dplyr)
load('../../01_DATA/ANGERS/LeafOptics.RData')
source('../Libraries/Lib_nrmse_func.R')
source('../Libraries/Lib_Analysis_Inversion.R')
################################################################################
# results data paths
################################################################################
pathdata1 <- "../../03_RESULTS/01_TEST_PAS_MAT/PAS="
pathdata3 <- "nm.RData"

################################################################################
# initialization
################################################################################
nrmse_dfLMA = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)
nrmse_dfEWT = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)

for (i in c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,175,200,225,250)) {
  pathdata2 <- as.character(i)
  pathdata01 <- paste(pathdata1,pathdata2,sep = "")
  pathdata02 <- paste(pathdata01,pathdata3, sep = "")
  load(pathdata02)
  
  NRMSELMA = as.numeric(get_performances_inversion(target = Biochemistry$LMA, 
                                                   estimate = parms_est_df$LMA, 
                                                   categories= TRUE)[3])
  NRMSEEWT = as.numeric(get_performances_inversion(target = Biochemistry$EWT, 
                                                   estimate = parms_est_df$EWT, 
                                                   categories= TRUE)[3])
  R2LMA = as.numeric(get_performances_inversion(target = Biochemistry$LMA, 
                                                   estimate = parms_est_df$LMA, 
                                                   categories= TRUE)[1])
  R2EWT = as.numeric(get_performances_inversion(target = Biochemistry$EWT, 
                                                   estimate = parms_est_df$EWT, 
                                                   categories= TRUE)[1])
  nrmse_dfLMA = rbind(nrmse_dfLMA,data.frame("var"= i, 
                                             "nrmse" = NRMSELMA, "r2" = R2LMA))
  nrmse_dfEWT = rbind(nrmse_dfEWT,data.frame("var"= i, 
                                             "nrmse" = NRMSEEWT, "r2" = R2EWT))
}

nrmse_dfLMA <- nrmse_dfLMA[-c(1),]
nrmse_dfEWT <- nrmse_dfEWT[-c(1),]

################################################################################
# Plots
################################################################################


plotLMA <- ggplot(nrmse_dfLMA, aes(x = var, y = nrmse)) +
  geom_line(aes(x = var, y = nrmse), colour = "black", size = 1) +
  geom_line(aes(x = var, y = (1-r2)*100), colour = "red", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (red) & 1-R² (black) (%)")+
  scale_x_continuous("STEP_LMA (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))

filename = file.path('../../03_RESULTS/01_TEST_PAS_MAT/NRMSE_LMA_standard_pas.png')
ggsave(filename,plot = plotLMA, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plotEWT <- ggplot(nrmse_dfEWT, aes(x = var, y = nrmse)) +
  geom_line(aes(x = var, y = nrmse), colour = "black", size = 1) +
  geom_line(aes(x = var, y =(1-r2)*100), colour = "red", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (red) & 1-R² (black) (%)")+
  scale_x_continuous("STEP_EWT (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))
filename = file.path('../../03_RESULTS/01_TEST_PAS_MAT/NRMSE_EWT_standard_pas.png')
ggsave(filename,plot = plotEWT, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

