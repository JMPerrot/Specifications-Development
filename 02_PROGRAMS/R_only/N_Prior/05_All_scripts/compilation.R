################################################################################
## This code summarizes the evolution of the performances of PROSPECT inversion
################################################################################
# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


## Libraries required ----------------------------------------------------------

library(tidyverse)
library(prospect)
library(data.table)
library(doFuture)
library(ggplot2)

source('../Libraries/Lib_Plots.R')
source('../Libraries/Lib_Analysis_Inversion.R')

## Scripts required ------------------------------------------------------------

# source("../../02_PROGRAMS/00_WriteDatabase/fun_01_Write_Database.r")
# source("../../02_PROGRAMS/01_Reference_PROSPECT_Inversion/fun_01_Baseline_estimations_REF#1.R")
# source("../../02_PROGRAMS/01_Reference_PROSPECT_Inversion/fun_01_Baseline_estimations_REF#2.R")
source("../../02_PROGRAMS/02_SpectralSampling/fun_01_SpectralSampling.R")
# source("../../02_PROGRAMS/03_Shift/fun_01_Shift.R")
#source("../../02_PROGRAMS/04_SFS/fun_01_SFS.R")

a<-SpecSampl()
