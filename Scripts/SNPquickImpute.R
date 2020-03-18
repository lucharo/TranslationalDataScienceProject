# SNP QUICK IMPUTE


library(impute)
library(tidyverse)

cluster = 1

if (cluster == 1){
  snp.original = readRDS("../FULLDATA/genetic_data_cvd_snps.rds")
  snp_info.original = readxl::read_xlsx("../SNP_info.xlsx")

  snp.original = cbind(ID = rownames(snp.original), snp.original)
  snp.original$ID = as.numeric(levels(snp.original$ID)[snp.original$ID])

  save_path = "../FULLDATA/preprocessed/"
} else {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

    snp.original = readRDS("../data/Genes_toy.rds")
    snp_info.original = readxl::read_xlsx("../SNP_info.xlsx")

    snp.original = cbind(ID = 1:nrow(snp.original), snp.original)
    save_path = "../data/preprocessed/"
}
getwd()

#################################################################
##                     Proccessing of SNPs                     ##
#################################################################

# removing any snps from the info file that are not included in the data provided by barbz
snps = colnames(snp.original)
snp.info = snp_info.original[snp_info.original$markername %in% snps, ]

saveRDS(snp.info, file = paste0(save_path,"snpInfo.rds"))


# recoding 00 as NA,  01 as 0, 02 as 1, and 03 as 2
snp = cbind(ID = snp.original$ID,
            as.data.frame(
              apply(snp.original[,-1], 2,
                    as.numeric)
            ))
# snp <- snp.original
# 
# for (i in colnames(snp)){
#   snp[,i] = as.numeric(snp[,i], multiple=TRUE)
# }



snp[snp==0] <- NA
snp[snp==1] <- 0
snp[snp==2] <- 1
snp[snp==3] <- 2

snp$ID = snp.original$ID

# SNP DATA A BIT HARDER TO IMPUTE BECAUSE VALUES ARE EITHER 0,1 OR 2, 
# SOME ENTRIES (PARTICUARLY COLUMN 77 HAS A LOT OF THE SAME VALUES, =2 FOR COLUMN 77)
# WHICH MAKES IT EASY TO DO IMPUTATION:
# -- mean(snp[,77]==2, na.rm = T) == 0.995998
# -- sum(colMeans(snp==0, na.rm = T)>0.7) = 0
# -- sum(colMeans(snp==1, na.rm = T)>0.95) = 0
# -- sum(colMeans(snp==2, na.rm = T)>0.95) = 9 # having these columns in results in a standard deviation of
# 0 which does not allow us to scale the SNP values, BUT MAYBE SNPs dont need to be rescaled
# -- sum(colMeans(snp==2, na.rm = T)>0.7) = 60
snp = snp[,!(colMeans(snp==2,na.rm= T)>0.9)]

# NOT REALLY SURE HOW TO IMPUTE SNPs, ASK BARBS

# RMSE.snps = t(sapply(1:10,function(x) kNNImputeOptimization(data.in = snp, seed = x)))
# boxplot(RMSE.snps)
# best.k.med.snps = which.min(apply(RMSE.snps, 2, median))
# best.k.mean.snps = which.min(colMeans(RMSE.snps))

#### FOR NOW, QUICK AND DIRTY SOLUTION
snp.imp = data.frame(round(impute.knn(as.matrix(snp[,-1]))$data,0))
snp.imp = cbind(ID = snp$ID, snp.imp)
snp.imp = as.data.frame(snp.imp, stringsAsFactors = FALSE)
stopifnot(all(snp.imp[,-1]%%1 == 0))
saveRDS(snp.imp, paste0(save_path,"snpImputed.rds"))

saveRDS(snp, paste0(save_path,"snpProcessed.rds"))
