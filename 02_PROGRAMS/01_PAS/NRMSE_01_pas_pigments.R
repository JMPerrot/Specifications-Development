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
pathdata1 <- "../../03_RESULTS/01_TEST_PAS/01_TEST_PAS_PIG/PAS="
extention <- "nm.RData"

################################################################################
# initialization
################################################################################
nrmse_dfCHL = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)
nrmse_dfCAR = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)


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
  R2CHL = as.numeric(get_performances_inversion(target = Biochemistry$CHLa+Biochemistry$CHLb, 
                                                   estimate = parms_est_df$CHL, 
                                                   categories= TRUE)[1])
  R2CAR = as.numeric(get_performances_inversion(target = Biochemistry$CAR, 
                                                   estimate = parms_est_df$CAR, 
                                                   categories= TRUE)[1])

  nrmse_dfCHL = rbind(nrmse_dfCHL,data.frame("var"= i, 
                                             "nrmse" = NRMSECHL, "r2" = R2CHL))
  nrmse_dfCAR = rbind(nrmse_dfCAR,data.frame("var"= i, 
                                             "nrmse" = NRMSECAR, "r2" = R2CAR))

}
nrmse_dfCHL <- nrmse_dfCHL[-c(1),]
nrmse_dfCAR <- nrmse_dfCAR[-c(1),]



plotCHL <- ggplot(nrmse_dfCHL, aes(x = var, y = nrmse)) +
  geom_line(aes(x = var, y = nrmse), colour = "black", size = 1) +
  geom_line(aes(x = var, y =(1-r2)*100), colour = "red", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (black -) & 1-R² (red .-) (%)")+
  scale_x_continuous("STEP_CHL (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))
filename = file.path('../../03_RESULTS/01_TEST_PAS/NRMSE_CHL_standard_pas.png')
ggsave(filename,plot = plotCHL, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)

plotCAR <- ggplot(nrmse_dfCAR, aes(x = var, y = nrmse)) +
  geom_line(aes(x = var, y = nrmse), colour = "black", size = 1) +
  geom_line(aes(x = var, y =(1-r2)*100), colour = "red", size = 1, linetype = "3313")+
  
  scale_y_continuous("NRMSE (black -) & 1-R² (red .-) (%)")+
  scale_x_continuous("STEP_CAR (nm)")+
  theme(
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red"))
filename = file.path('../../03_RESULTS/01_TEST_PAS/NRMSE_CAR_standard_pas.png')
ggsave(filename,plot = plotCAR, device = "png", path = NULL,
       scale = 1, width = 20, height = 13, units = "cm",
       dpi = 600)