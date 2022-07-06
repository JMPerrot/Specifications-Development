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

source('../../Libraries/Lib_Plots.R')
source('../../Libraries/Lib_Analysis_Inversion.R')

## Scripts required ------------------------------------------------------------

# source("../00_WriteDatabase/01_Write_Database.r")
source("../01_Reference_PROSPECT_Inversion/01_Baseline_estimations_REF#1.R")
source("../01_Reference_PROSPECT_Inversion/01_Baseline_estimations_REF#2.R")

source("../02_SpectralSampling/main_01_SpectralSampling.R")
source("../02_SpectralSamplingmain_02_plot_SpectralSampling.R")

source("../03_Shift/main_01_Shift.R")
source("../03_Shift/main_02_plot_SpectralShifting.R")

source("../04_SFS/main_01_SFS.R")


