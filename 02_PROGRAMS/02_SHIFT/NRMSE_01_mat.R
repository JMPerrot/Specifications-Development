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
shift_max = 200

pathdata1 <- paste("../../03_RESULTS/02_TEST_SHIFT/pas=",shift_max,"nm/SHIFT=", sep = "")
extention <- "nm.RData"

################################################################################
# initialization
################################################################################
nrmse_dfEWT = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)
nrmse_dfLMA = data.frame("var"=NA, "nrmse" = NA, "r2"= NA)



for (i in c(1:shift_max)) {
  pas <- as.character(i)
  pathdata01 <- paste(pathdata1,pas,sep = "")
  pathdata02 <- paste(pathdata01,extention, sep = "")
  load(pathdata02)
  
  NRMSEEWT = as.numeric(get_performances_inversion(target = Biochemistry$EWT, 
                                                   estimate = parms_est_df$EWT, 
                                                   categories= TRUE)[3])
  NRMSELMA = as.numeric(get_performances_inversion(target = Biochemistry$LMA, 
                                                   estimate = parms_est_df$LMA, 
                                                   categories= TRUE)[3])
  R2EWT = as.numeric(get_performances_inversion(target = Biochemistry$EWTa, 
                                                estimate = parms_est_df$EWT, 
                                                categories= TRUE)[1])
  R2LMA = as.numeric(get_performances_inversion(target = Biochemistry$LMA, 
                                                estimate = parms_est_df$LMA, 
                                                categories= TRUE)[1])
  
  nrmse_dfEWT = rbind(nrmse_dfEWT,data.frame("var"= i, 
                                             "nrmse" = NRMSEEWT, "r2" = R2EWT))
  nrmse_dfLMA = rbind(nrmse_dfLMA,data.frame("var"= i, 
                                             "nrmse" = NRMSELMA, "r2" = R2LMA))
  
}
nrmse_dfEWT <- nrmse_dfEWT[-c(1),]
nrmse_dfLMA <- nrmse_dfLMA[-c(1),]
nrmse_mat<-data.frame("EWT" <- nrmse_dfEWT, "LMA" <- nrmse_dfLMA)
save(nrmse_mat, file = paste("../../03_RESULTS/02_TEST_SHIFT/pas=",shift_max, "nm/nrmse_mat.RData", sep = ""))


pEWT<-ggplot(nrmse_dfEWT, aes(x = var, y = nrmse))+
  geom_line(aes(x = var, y = nrmse), size =1, color = "blue")+
  xlab("translation (nm)")+
  ylab("NRMSE (%)")+    
  annotate("text", x = 10, y = min(nrmse_dfEWT$nrmse), label = paste("pas_",as.character(shift_max),"nm", sep = ""),
           hjust = 0, parse = TRUE)

ggsave(filename = paste("../../03_RESULTS/02_TEST_SHIFT/EWT_shift=",shift_max,"nm.png",sep = ""),plot = pEWT)

pLMA<-ggplot(nrmse_dfLMA, aes(x = var, y = nrmse))+
  geom_line(aes(x = var, y = nrmse), size =1, color = "red")+
  xlab("translation (nm)")+
  ylab("NRMSE (%)")+    
  annotate("text", x = 10, y = min(nrmse_dfLMA$nrmse), label = paste("pas_",as.character(shift_max),"nm", sep = ""),
           hjust = 0, parse = TRUE)

ggsave(filename = paste("../../03_RESULTS/02_TEST_SHIFT/LMA_shift=",shift_max,"nm.png",sep = ""),plot = pLMA)
