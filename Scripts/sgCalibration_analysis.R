
library(ggplot2)
library(dplyr)
library(ggpubr)

platform = Sys.info()['sysname']
if (platform == 'Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/PLS"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

bio = readRDS(paste0(data_folder,"bioImputedKNN.rds"))[,-1]

read_data = paste0(save_plots,"ArrayJob/")


sg_calib = readRDS(paste0(read_data, "sgCalibration_1"))
quietly(sapply(2:100, FUN =  function(x) {
  assign("sg_calib", rbind(sg_calib, readRDS(paste0(read_data, "sgCalibration_",x,".rds"))),
         envir = .GlobalEnv)}))
lasso.stab = data.frame(lasso.stab)
colnames(lasso.stab) = colnames(bio)

