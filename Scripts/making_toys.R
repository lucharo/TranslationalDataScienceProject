# Making toy dataset(s) from UK biobank

library(clusteringdatasets)
setwd("/rds/general/user/lc5415/home/hda_tds_ukbiobank")
b = readRDS("Biomarkers.rds")
c = readRDS("Covariates.rds")

L = 1000
b.toy = b[sample(c(1:nrow(b)), size = L),]
c.toy = c[sample(c(1:nrow(c)), size = L),]

saveRDS(b.toy, file = "bio_toy.rds")
saveRDS(c.toy, file = "cov_toy.rds")
