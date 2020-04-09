# KNN imputation of bio
library(impute)
library(tidyverse)
library(ggplot2)

cluster = 1

if (cluster == 1){
 save_data = data_folder = "../FULLDATA/preprocessed/"
 save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}


ifelse(!dir.exists(file.path(save_plots, "knnOpti/")),
       dir.create(file.path(save_plots, "knnOpti/")), FALSE)
save_plots = paste0(save_plots,"knnOpti/")

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))
# log.bio = log10(bio)
# saveRDS(log.bio, paste0(save_data,"LOGbioProcessed.rds"))
log.bio = readRDS(paste0(data_folder,"LOGbioProcessed.rds"))

print("Starting Imputation #1...")
# Impute with KNN -- using caret
t0 = Sys.time()

source("kNNImputeOptimization.R",print.eval = T)

############################################################################
############################################################################
###                                                                      ###
###                          BIO KNN IMPUTATION                          ###
###                                                                      ###
############################################################################
############################################################################

#### take in command args
args = commandArgs(trailingOnly = TRUE)

CV = as.numeric(args[1])

print(CV)

results = t(sapply(CV,
                      function(x) kNNImputeOptimization(bio, seed = x,
                                                  perParam = T, scaled = T,
                                                  plot = F))) 
                      
RMSE = results[[1]]

RMSE = t(RMSE)
                
ifelse(!dir.exists(file.path(save_plots, "ArrayJob/")),
       dir.create(file.path(save_plots, "ArrayJob/")), FALSE)
save_plots = paste0(save_plots,"ArrayJob/")

saveRDS(RMSE, paste0(save_plots, "lassoStab",as.character(CV),".rds"))
                   

print(Sys.time() - t0) # takes about 1 minute
print("Imputation number 1 finished.")