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

source("kNNImputeOptimization.R",print.eval = T)

############################################################################
############################################################################
###                                                                      ###
###                          BIO KNN IMPUTATION                          ###
###                                                                      ###
############################################################################
############################################################################

CV = 10
results = t(sapply(1:CV,
                function(x) kNNImputeOptimization(bio, seed = x,
                                                  perParam = T, scaled = T,
                                                  plot = x==CV)))
RMSE = results[[1]]
for (i in 2:CV){
  if (i != CV){
    RMSE = cbind(RMSE, results[[i]])
  } else{
    RMSE = cbind(RMSE, results[[CV]][[1]])
  }
}

RMSE = t(RMSE)
factplot = results[[CV]][[2]]
bplot = results[[CV]][[3]]
scatplot = results[[CV]][[4]]

pdf(paste0(save_plots,"bioFacetNA.pdf"))
factplot
dev.off()

pdf(paste0(save_plots,"bioBarNA.pdf"))
bplot
dev.off()

pdf(paste0(save_plots,"bioScatNA.pdf"))
scatplot
dev.off()

pdf(paste0(save_plots,"bioFacetNA.pdf"))
boxplot(RMSE)
title(main = paste0("log bio RMSE per choice of k,\n min: ",
                    as.character(min(colMeans(RMSE))),
      "\n k: ", as.character(which.min(colMeans(RMSE)))))
dev.off()

boxplot(RMSE)
best.k.med = which.min(apply(RMSE, 2, median))
best.k.mean = which.min(colMeans(RMSE))

# ?scaling: works weird man
mean.cols = as.vector(colMeans(bio, na.rm = T))
sd.cols = as.vector(apply(bio,2, function(col) sd(col, na.rm = T)))
bio.scaled = sweep(sweep(bio,MARGIN = 2,mean.cols,'-'),2,sd.cols,"/")
bio.scaled = as.matrix(bio.scaled)
bio.scaled.imp = impute.knn(bio.scaled, k = best.k.mean)$data

#descale 0r rescale, however you wanna call it
bio.imp = sweep(sweep(bio.scaled.imp,MARGIN = 2,sd.cols,'*'),2,mean.cols,"+")
# this weird
# matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T)*c(2,3,2)
# this works
#sweep(matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T), 2, 
# c(2,3,2), "*")
print(Sys.time() - t0) # takes about 1 minute


saveRDS(bio.imp, file = paste0(save_data,"bioImputed.rds"))

############################################################################
############################################################################
###                                                                      ###
###                        LOG BIO KNN IMPUTATION                        ###
###                                                                      ###
############################################################################
############################################################################

results = t(sapply(1:CV,function(x) kNNImputeOptimization(log10(bio), seed = x,
                                                      perParam = T, scaled = T,
                                                      plot = x == CV)))
RMSE = results[[1]]
for (i in 2:CV){
  if (i != CV){
    RMSE = cbind(RMSE, results[[i]])
  } else{
    RMSE = cbind(RMSE, results[[CV]][[1]])
  }
}
RMSE = t(RMSE)
factplot = results[[CV]][[2]]
bplot = results[[CV]][[3]]
scatplot = results[[CV]][[4]]

pdf(paste0(save_plots,"logbioFacetNA.pdf"))
factplot
dev.off()

pdf(paste0(save_plots,"logbioBarNA.pdf"))
bplot
dev.off()

pdf(paste0(save_plots,"logbioScatNA.pdf"))
scatplot
dev.off()

pdf(paste0(save_plots,"logbioFacetNA.pdf"))
boxplot(RMSE)
title(main = paste0("log bio RMSE per choice of k,\n min: ",
                    as.character(min(colMeans(RMSE))),
                    "\n k: ", as.character(which.min(colMeans(RMSE)))))
dev.off()

best.k.med = which.min(apply(RMSE, 2, median))
best.k.mean = which.min(colMeans(RMSE))

# ?scaling: works weird man
mean.cols = as.vector(colMeans(bio, na.rm = T))
sd.cols = as.vector(apply(bio,2, function(col) sd(col, na.rm = T)))
bio.scaled = sweep(sweep(bio,MARGIN = 2,mean.cols,'-'),2,sd.cols,"/")
bio.scaled = as.matrix(bio.scaled)
bio.scaled.imp = impute.knn(bio.scaled, k = best.k.mean)$data

#descale 0r rescale, however you wanna call it
bio.imp = sweep(sweep(bio.scaled.imp,MARGIN = 2,sd.cols,'*'),2,mean.cols,"+")
# this weird
# matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T)*c(2,3,2)
# this works
#sweep(matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T), 2, 
# c(2,3,2), "*")
print(Sys.time() - t0) # takes about 1 minute


saveRDS(bio.imp, file = paste0(save_data,"LOGbioImputed.rds"))