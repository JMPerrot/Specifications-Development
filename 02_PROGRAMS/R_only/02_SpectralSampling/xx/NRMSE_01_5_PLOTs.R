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

################################################################################
# Data loading 
################################################################################

load("../../03_RESULTS/01_TEST_PAS/nrmse_mat.RData")
load("../../03_RESULTS/01_TEST_PAS/nrmse_pig.RData")

################################################################################
# Formatting
################################################################################

nrmse_dfCHL<-data.frame("wl" = nrmse_pig[,1], "nrmse" = nrmse_pig[,2], "r2" = nrmse_pig[,3])
nrmse_dfCAR<-data.frame("wl" = nrmse_pig[,4], "nrmse" = nrmse_pig[,5], "r2" = nrmse_pig[,6])

nrmse_dfEWT<-data.frame("wl" = nrmse_mat[,1], "nrmse" = nrmse_mat[,2], "r2" = nrmse_mat[,3])
nrmse_dfLMA<-data.frame("wl" = nrmse_mat[,4], "nrmse" = nrmse_mat[,5], "r2" = nrmse_mat[,6])
################################################################################
# Plots
################################################################################


plotLMA <- ggplot(nrmse_dfLMA, aes(x = wl, y = nrmse)) +
  geom_line(aes(x = wl, y = nrmse), colour = "#000000", size = 1) +
  geom_line(aes(x = wl, y = (1-r2)*100), colour = "#999999", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (dark -) & 1-R² (light .-) (%)")+
  scale_x_continuous("STEP_LMA (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"))


filename = file.path('../../03_RESULTS/01_TEST_PAS/NRMSE_LMA_standard_pas.png')
ggsave(filename,plot = plotLMA, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plotEWT <- ggplot(nrmse_dfEWT, aes(x = wl, y = nrmse)) +
  geom_line(aes(x = wl, y = nrmse), colour = "#3300CC", size = 1) +
  geom_line(aes(x = wl, y =(1-r2)*100), colour = "#66CCFF", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (dark -) & 1-R² (light .-) (%)")+
  scale_x_continuous("STEP_EWT (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"))

filename = file.path('../../03_RESULTS/01_TEST_PAS/NRMSE_EWT_standard_pas.png')
ggsave(filename,plot = plotEWT, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plotCHL <- ggplot(nrmse_dfCHL, aes(x = wl, y = nrmse)) +
  geom_line(aes(x = wl, y = nrmse), colour = "#009900", size = 1) +
  geom_line(aes(x = wl, y =(1-r2)*100), colour = "#00FF33", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (dark -) & 1-R² (light .-) (%)")+
  scale_x_continuous("STEP_CHL (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"))

filename = file.path('../../03_RESULTS/01_TEST_PAS/NRMSE_CHL_standard_pas.png')
ggsave(filename,plot = plotCHL, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plotCAR <- ggplot(nrmse_dfCAR, aes(x = wl, y = nrmse)) +
  geom_line(aes(x = wl, y = nrmse), colour = "#FF6600", size = 1) +
  geom_line(aes(x = wl, y =(1-r2)*100), colour = "#FFCC00", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (dark -) & 1-R² (light .-) (%)")+
  scale_x_continuous("STEP_CAR (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"))

filename = file.path('../../03_RESULTS/01_TEST_PAS/NRMSE_CAR_standard_pas.png')
ggsave(filename,plot = plotCAR, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plot_nrmse1 <- ggarrange(plotCHL,plotCAR,plotEWT,plotLMA,
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
plot_nrmse<- annotate_figure(plot_nrmse1, 
                           bottom = text_grob("NRMSE et (1-R²) en fonction du pas", 
                                              color = "black", 
                                              face = "bold", 
                                              size = 14))

ggsave("../../03_RESULTS/01_TEST_PAS/stats_comp.png", plot_nrmse,device = "png")
