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
load('../../01_DATA/ANGERS/LeafOptics.RData')
source('../Libraries/Lib_nrmse_func.R')
source('../Libraries/Lib_Analysis_Inversion.R')
################################################################################
# results data paths
################################################################################
pathdata1 <- "../../03_RESULTS/01_TEST_PAS_PIG/PAS="
extention <- "nm.RData"

################################################################################
# initialization
################################################################################
nrmse_dfCHL = data.frame("var"=NA, "nrmse" = NA)
nrmse_dfCAR = data.frame("var"=NA, "nrmse" = NA)
nrmse_dfLMA = data.frame("var"=NA, "nrmse" = NA)
nrmse_dfEWT = data.frame("var"=NA, "nrmse" = NA)

for (i in c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,175,200,225,250)) {
  pas <- as.character(i)
  pathdata01 <- paste(pathdata1,pas,sep = "")
  pathdata02 <- paste(pathdata01,extention, sep = "")
  load(pathdata02)
  
  NRMSECHL = as.numeric(get_performances_inversion(target = Biochemistry$CHLa+Biochemistry$CHLb, 
                                                   estimate = parms_est_df$CHL, 
                                                   categories= TRUE)[3])
  NRMSECAR = as.numeric(get_performances_inversion(target = Biochemistry$CAR, 
                                                   estimate = parms_est_df$CAR, 
                                                   categories= TRUE)[3])
  NRMSELMA = as.numeric(get_performances_inversion(target = Biochemistry$LMA, 
                                                   estimate = parms_est_df$LMA, 
                                                   categories= TRUE)[3])
  NRMSEEWT = as.numeric(get_performances_inversion(target = Biochemistry$EWT, 
                                                   estimate = parms_est_df$EWT, 
                                                   categories= TRUE)[3])
  nrmse_dfCHL = rbind(nrmse_dfCHL,data.frame("var"= i, 
                                             "nrmse" = NRMSECHL))
  nrmse_dfCAR = rbind(nrmse_dfCAR,data.frame("var"= i, 
                                             "nrmse" = NRMSECAR))
  nrmse_dfLMA = rbind(nrmse_dfLMA,data.frame("var"= i, 
                                             "nrmse" = NRMSELMA))
  nrmse_dfEWT = rbind(nrmse_dfEWT,data.frame("var"= i, 
                                             "nrmse" = NRMSEEWT))
}
nrmse_dfCHL <- nrmse_dfCHL[-c(1),]
nrmse_dfCAR <- nrmse_dfCAR[-c(1),]
nrmse_dfLMA <- nrmse_dfLMA[-c(1),]
nrmse_dfEWT <- nrmse_dfEWT[-c(1),]


plotCHL <- ggplot(nrmse_dfCHL, aes(x = var, y = nrmse)) +
  geom_point(aes(x = var, y = nrmse), colour = "#66CC00", size = 2) +
  geom_line(aes(x = var, y = nrmse), colour = "#66CC00", size = 1) +
  scale_size_manual(values = c(5, 5)) +
  labs(x = 'CHL_PAS (nm)',y = 'NRMSE') +
  ylim(0,max(nrmse_dfCHL$nrmse)) +
  xlim(0,max(nrmse_dfCHL$var))
filename = file.path('../../03_RESULTS/01_TEST_PAS_PIG/NRMSE_CHL_standard_pas.png')
ggsave(filename,plot = plotCHL, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plotCAR <- ggplot(nrmse_dfCAR, aes(x = var, y = nrmse)) +
  geom_point(aes(x = var, y = nrmse), colour = "orange", size = 2) +
  geom_line(aes(x = var, y = nrmse), colour = "orange", size = 1) +
  scale_size_manual(values = c(5, 5)) +
  labs(x = 'CAR_PAS (nm)',y = 'NRMSE') +
  ylim(0,max(nrmse_dfCAR$nrmse)) +
  xlim(0,max(nrmse_dfCAR$var))
filename = file.path('../../03_RESULTS/01_TEST_PAS_PIG/NRMSE_CAR_standard_pas.png')
ggsave(filename,plot = plotCAR, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

# plotLMA <- ggplot(nrmse_dfLMA, aes(x = var, y = nrmse)) +
#   geom_point(aes(x = var, y = nrmse), colour = "red", size = 2) +
#   geom_line(aes(x = var, y = nrmse), colour = "red", size = 1) +
#   scale_size_manual(values = c(5, 5)) +
#   labs(x = 'LMA_PAS (nm)',y = 'NRMSE') +
#   ylim(0,max(nrmse_dfLMA$nrmse)) +
#   xlim(0,max(nrmse_dfLMA$var))
# filename = file.path('../../03_RESULTS/02_TEST_PAS/NRMSE_LMA_standard_pas.png')
# ggsave(filename,plot = plotLMA, device = "png", path = NULL,
#        scale = 1, width = 20, height = 13, units = "cm",
#        dpi = 600)
# 
# plotEWT <- ggplot(nrmse_dfEWT, aes(x = var, y = nrmse)) +
#   geom_point(aes(x = var, y = nrmse), colour = "red", size = 2) +
#   geom_line(aes(x = var, y = nrmse), colour = "red", size = 1) +
#   scale_size_manual(values = c(5, 5)) +
#   labs(x = 'EWT_PAS (nm)',y = 'NRMSE') +
#   ylim(0,max(nrmse_dfEWT$nrmse)) +
#   xlim(0,max(nrmse_dfEWT$var))
# filename = file.path('../../03_RESULTS/02_TEST_PAS/NRMSE_EWT_standard_pas.png')
# ggsave(filename,plot = plotEWT, device = "png", path = NULL,
#        scale = 1, width = 20, height = 13, units = "cm",
#        dpi = 600)
# 
# 
# 
