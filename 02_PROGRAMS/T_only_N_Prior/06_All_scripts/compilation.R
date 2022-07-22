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

# source("../01_Reference_PROSPECT_Inversion/01_Baseline_estimations_REF#1_T_N.R")
# source("../01_Reference_PROSPECT_Inversion/01_Baseline_estimations_REF#2_T_N.R")

source("../02_SpectralSampling/main_01_SpectralSampling_T_N.R")
source("../02_SpectralSampling/main_02_plot_SpectralSampling_T_N.R")

source("../03_Shift/main_01_Shift_T_N.R")
source("../03_Shift/main_02_plot_SpecShifting_T_N.R")

source("../04_SFS/main_01_SFS_T_N.R")


