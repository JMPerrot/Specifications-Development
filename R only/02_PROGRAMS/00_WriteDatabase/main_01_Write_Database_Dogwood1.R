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
PathData <- '../../01_DATA'
dbName <- 'Dogwood1'
PathLOPdb <- file.path(PathData,dbName)
dir.create(path = PathLOPdb,showWarnings = F,recursive = T)

################################################################################
# download leaf optical properties and chemistry from online repository
################################################################################
# repository where data are stored
path_DataBioch<-file.path(PathLOPdb,"DataBioch.txt")
bioch<-read.table(path_DataBioch, header = FALSE, sep = "", dec = ".")
Biochemistry<-data.frame("CHL" = bioch$V1, "CAR" = bioch$V2, "EWT" = bioch$V3, "LMA" = bioch$V4)
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
