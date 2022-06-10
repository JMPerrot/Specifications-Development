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
dbName <- 'ANGERS'
PathLOPdb <- file.path(PathData,dbName)
dir.create(path = PathLOPdb,showWarnings = F,recursive = T)

################################################################################
# download leaf optical properties and chemistry from online repository
################################################################################
# repository where data are stored
gitlab_Rep <- 'https://gitlab.com/jbferet/myshareddata/raw/master/LOP'
# files available
fileName <- list(Refl = 'ReflectanceData.txt',
                 Tran = 'TransmittanceData.txt',
                 Bioch = 'DataBioch.txt')
Reflectance <- data.table::fread(file.path(gitlab_Rep,dbName,fileName$Refl))
Transmittance <- data.table::fread(file.path(gitlab_Rep,dbName,fileName$Tran))
Biochemistry <- data.table::fread(file.path(gitlab_Rep,dbName,fileName$Bioch))
Biochemistry$CHL <- Biochemistry$CHLa + Biochemistry$CHLb

List2Save <- list('Reflectance'=Reflectance,
                  'Transmittance'=Transmittance,
                  'Biochemistry'=Biochemistry)
path_save <- file.path(PathLOPdb,'LeafOptics.RData')
save(Reflectance,Transmittance,Biochemistry,file = path_save)
