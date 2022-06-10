################################################################################
## This program aims downloading a dataset including leaf optical properties  ##
## from a gitlab repository, and saving it on local drive for next processings##
################################################################################
# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

################################################################################
# Define input output directories
################################################################################
PathData <- '../../../01_DATA'



################################################################################
# download leaf optical properties and chemistry from online repository
################################################################################
# repository where data are stored
data_base <- c('Dogwood1', 'Maple', 'Partenocissus')

for (name in data_base) {
  dbName <- name
  PathLOPdb <- file.path(PathData,dbName)
  # dir.create(path = PathLOPdb,showWarnings = F,recursive = T)
  path_DataBioch<-file.path(PathLOPdb,"DataBioch.txt")
  
  bioch<-read.table(path_DataBioch, header = FALSE, sep = "", dec = ".")
  Biochemistry<-data.frame("CHL" = bioch$V1+bioch$V2, "CAR" = bioch$V3, "ANT" = bioch$V4)
  
  path_Reflectance<-file.path(PathLOPdb,"ReflectanceData.txt")
  Reflectance<-read.table(path_Reflectance, header = FALSE, sep = "", dec = ".")
 
   path_Transmittance<-file.path(PathLOPdb,"TransmittanceData.txt")
  Transmittance<-read.table(path_Transmittance, header = FALSE, sep = "", dec = ".")
  # files available
  List2Save <- list('Reflectance'=Reflectance,
                    'Transmittance'=Transmittance,
                    'Biochemistry'=Biochemistry)
  path_save <- file.path(PathLOPdb,'LeafOptics.RData')
  save(Reflectance,Transmittance,Biochemistry,file = path_save)
}

################################################################################
#
################################################################################


Path_Data1 <- file.path(PathData,'Dogwood1','LeafOptics.RData')
load(Path_Data1)
Dogwood1_Biochemistry <- data.frame("Biochemistry" = Biochemistry)
Dogwood1_R <- data.frame("Reflectance" = Reflectance)
Dogwood1_T <- data.frame("Transmittance"  = Transmittance)
R_Dogwood1 <- Dogwood1_R[c(0:345),]
T_Dogwood1 <- Dogwood1_T[c(0:345),]

Path_Data2 <- file.path(PathData,'Maple','LeafOptics.RData')
load(Path_Data2)
Maple_Biochemistry <- data.frame("Biochemistry" = Biochemistry)
Maple_R <- data.frame("Reflectance" = Reflectance[,-1])
Maple_T <-data.frame("Transmittance"  = Transmittance[,-1])
R_Maple <- Maple_R[c(37:381),]
T_Maple <- Maple_T[c(37:381),]

Path_Data3 <- file.path(PathData,'Partenocissus','LeafOptics.RData')
load(Path_Data3)
Partenocissus_Biochemistry <- data.frame("Biochemistry" = Biochemistry)
Partenocissus_R <- data.frame("Reflectance" = Reflectance[,-1])
Partenocissus_T <- data.frame("Transmittance"  = Transmittance[,-1])
R_Partenocissus<-Partenocissus_R[c(37:381),]
T_Partenocissus<-Partenocissus_T[c(37:381),]

Reflectance <- unname(cbind(R_Dogwood1,R_Partenocissus,R_Maple))
Transmittance <- unname(cbind(T_Dogwood1,T_Partenocissus,T_Maple))
Biochemistry <- rbind(Partenocissus_Biochemistry,Maple_Biochemistry,Dogwood1_Biochemistry)

Savedatapath <- '../../../01_DATA/Test_data_LeafOptics.RData'
save(Reflectance,Transmittance,Biochemistry,file = Savedatapath)
load(Savedatapath)
