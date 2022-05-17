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
source('../Libraries/Lib_Analysis_Inversion.R')
################################################################################
# results data paths
################################################################################
shift_max = 150

pathdata1 <- paste("../../03_RESULTS/02_TEST_SHIFT/pas=",shift_max,"nm/SHIFT=", sep = "")
extention <- "nm.RData"

################################################################################
# initialization
################################################################################
nrmse_dfCHL = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)
nrmse_dfCAR = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)



for (i in c(1:shift_max)) {
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
nrmse_pig<-data.frame("CHL" <- nrmse_dfCHL, "CAR" <- nrmse_dfCAR)
save(nrmse_pig, file = paste("../../03_RESULTS/02_TEST_SHIFT/pas=",shift_max, "nm/nrmse_pig.RData", sep = ""))


pCHL<-ggplot(nrmse_dfCHL, aes(x = var, y = nrmse))+
  geom_line(aes(x = var, y = nrmse), size =1, color = "#009900")+
  xlab("translation (nm)")+
  ylab("NRMSE (%)")+    
  annotate("text", x = 10, y = min(nrmse_dfCHL$nrmse), label = paste("pas_",as.character(shift_max),"nm", sep = ""),
           hjust = 0, parse = TRUE)

ggsave(filename = paste("../../03_RESULTS/02_TEST_SHIFT/CHL_shift=",shift_max,"nm.png",sep = ""),plot = pCHL)

pCAR<-ggplot(nrmse_dfCAR, aes(x = var, y = nrmse))+
  geom_line(aes(x = var, y = nrmse), size =1, color = "orange")+
  xlab("translation (nm)")+
  ylab("NRMSE (%)")+    
  annotate("text", x = 10, y = min(nrmse_dfCAR$nrmse), label = paste("pas_",as.character(shift_max),"nm", sep = ""),
           hjust = 0, parse = TRUE)

ggsave(filename = paste("../../03_RESULTS/02_TEST_SHIFT/CAR_shift=",shift_max,"nm.png",sep = ""),plot = pCAR)
