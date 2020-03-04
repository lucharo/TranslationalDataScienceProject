# KNN imputation of bio
library(impute)
cluster = 0
data_folder = "../data/preprocessed/"
if (cluster == 1){
 data_folder = "../FULLDATA/preprocessed/"
}

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))

# Impute with KNN -- using caret
t0 = Sys.time()

source("kNNImputeOptimization.R",print.eval = T)

RMSE = t(sapply(1:5,
                function(x) kNNImputeOptimization(bio, seed = x,
                                                  perParam = T, scaled = T,
                                                  plot = x==5)))
boxplot(RMSE)
best.k.med = which.min(apply(RMSE, 2, median))
best.k.mean = which.min(colMeans(RMSE))

RMSE = t(sapply(1:5,function(x) kNNImputeOptimization(log10(bio), seed = x)))
boxplot(RMSE)
best.k.med = which.min(apply(RMSE, 2, median))
best.k.mean = which.min(colMeans(RMSE))

# ?scaling: works weird man
mean.cols = as.vector(colMeans(bio, na.rm = T))
sd.cols = as.vector(apply(bio,2, function(col) sd(col, na.rm = T)))
bio.scaled = sweep(sweep(bio,MARGIN = 2,mean.cols,'-'),2,sd.cols,"/")
bio.scaled = as.matrix(bio.scaled)
bio.scaled.imp = impute.knn(bio.scaled)$data

#descale 0r rescale, however you wanna call it
bio.imp = sweep(sweep(bio.scaled.imp,MARGIN = 2,sd.cols,'*'),2,mean.cols,"+")
# this weird
# matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T)*c(2,3,2)
# this works
#sweep(matrix(c(1,2,3,5,3,4), nrow = 2, byrow = T), 2, 
# c(2,3,2), "*")
print(Sys.time() - t0) # takes about 1 minute


saveRDS(bio.imp, file = paste0(save_path,"bioImputed.rds"))