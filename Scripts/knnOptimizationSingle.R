print(.libPaths())

# KNN imputation of bio
library(impute)
if (!require("tidyverse")) install.packages("tidyverse",  repos='http://cran.us.r-project.org')
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
bio = bio %>% select(-ID)
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

results = kNNImputeOptimization(bio, seed = CV,
                                perParam = T, scaled = T,
                                plot = F) 
                   
print(typeof(results))
print(results)
                      
ifelse(!dir.exists(file.path(save_plots, "ArrayJob/")),
       dir.create(file.path(save_plots, "ArrayJob/")), FALSE)
save_plots = paste0(save_plots,"ArrayJob/")

saveRDS(results, paste0(save_plots, "knnOpti",as.character(CV),".rds"))
                   

print(Sys.time() - t0) # takes about 1 minute