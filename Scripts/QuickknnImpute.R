# KNN imputation of bio
library(impute)
library(tidyverse)
library(ggplot2)

cluster = 0

if (cluster == 1){
 save_data = data_folder = "../FULLDATA/preprocessed/"
 save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))
log.bio = as.data.frame(cbind(ID = bio$ID, log10(bio[,-1])),
                        stringsASFactors = F)

saveRDS(log.bio, paste0(save_data,"LOGbioProcessed.rds"))

## remove IDs, add later
ids = bio$ID
bio = bio[,-1]
log.bio = log.bio[,-1]

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

bio.imp = as.data.frame(cbind(ID = ids, bio.imp), stringsAsFactors = F)
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
mean.cols = as.vector(colMeans(log.bio, na.rm = T))
sd.cols = as.vector(apply(log.bio,2, function(col) sd(col, na.rm = T)))
log.bio.scaled = sweep(sweep(log.bio,MARGIN = 2,mean.cols,'-'),2,sd.cols,"/")
log.bio.scaled = as.matrix(log.bio.scaled)
log.bio.scaled.imp = impute.knn(log.bio.scaled, k = 15)$data

#descale 0r rescale, however you wanna call it
log.bio.imp = sweep(sweep(log.bio.scaled.imp,MARGIN = 2,sd.cols,'*'),2,mean.cols,"+")
# this weird
# matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T)*c(2,3,2)
# this works
#sweep(matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T), 2, 
# c(2,3,2), "*")
print(Sys.time() - t0) # takes about 1 minute

log.bio.imp = as.data.frame(cbind(ID = ids, log.bio.imp), stringsAsFactors = F)
saveRDS(log.bio.imp, file = paste0(save_data,"LOGbioImputedKNN.rds"))