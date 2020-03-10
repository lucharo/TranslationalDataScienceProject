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

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))
log.bio = log10(bio)
saveRDS(log.bio, paste0(save_data,"LOGbioProcessed.rds"))

# Impute with KNN -- using caret
t0 = Sys.time()

# ?scaling: works weird man
mean.cols = as.vector(colMeans(bio, na.rm = T))
sd.cols = as.vector(apply(bio,2, function(col) sd(col, na.rm = T)))
bio.scaled = sweep(sweep(bio,MARGIN = 2,mean.cols,'-'),2,sd.cols,"/")
bio.scaled = as.matrix(bio.scaled)
bio.scaled.imp = impute.knn(bio.scaled, k = 15)$data

#descale 0r rescale, however you wanna call it
bio.imp = sweep(sweep(bio.scaled.imp,MARGIN = 2,sd.cols,'*'),2,mean.cols,"+")
# this weird
# matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T)*c(2,3,2)
# this works
#sweep(matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T), 2, 
# c(2,3,2), "*")
                                      
print(Sys.time() - t0) # takes about 1 minute


saveRDS(bio.imp, file = paste0(save_data,"bioImputedKNN.rds"))

############################################################################
############################################################################
###                                                                      ###
###                        LOG BIO KNN IMPUTATION                        ###
###                                                                      ###
############################################################################
############################################################################

t0 = Sys.time()
# ?scaling: works weird man
mean.cols = as.vector(colMeans(bio, na.rm = T))
sd.cols = as.vector(apply(bio,2, function(col) sd(col, na.rm = T)))
bio.scaled = sweep(sweep(bio,MARGIN = 2,mean.cols,'-'),2,sd.cols,"/")
bio.scaled = as.matrix(bio.scaled)
bio.scaled.imp = impute.knn(bio.scaled, k = 15)$data

#descale 0r rescale, however you wanna call it
bio.imp = sweep(sweep(bio.scaled.imp,MARGIN = 2,sd.cols,'*'),2,mean.cols,"+")
# this weird
# matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T)*c(2,3,2)
# this works
#sweep(matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T), 2, 
# c(2,3,2), "*")
print(Sys.time() - t0) # takes about 1 minute


saveRDS(bio.imp, file = paste0(save_data,"LOGbioImputedKNN.rds"))