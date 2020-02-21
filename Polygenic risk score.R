# Aim of this script is to calculate a polygenic risk score for each person

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

snp.info <- readRDS("data/preprocessed/snpInfo.rds")
snp <- readRDS("data/preprocessed/snpProcessed.rds")


##################################################################
##                        Computing PRS                         ##
##################################################################

# Check that the snps are in the same order in snp and snp.info
stopifnot(all(colnames(snp) == snp.info$markername))

### Function for computing PRS per person
betas <- snp.info$beta

#this for loop calculates beta*copy number for each snp  
for (i in length(snp)){
  PRStest = snp[i,]*betas
}
#putting this in a function
snps.betas <- function() {
  for (i in length(snp)){
    PRStest = snp[i,]*betas
  }
}

#now including the weighted sum for each person (NAs are still present in the snp data so have included na.rm=TRUE for now, can remove if we do imputation)
compute_PRS <- function(data=snp){
  for (i in length(snp)){
    PRS_int = snp[i,]*betas
  }
  PRS = rowSums(PRS_int, na.rm=TRUE)
}

PRSt = snp[1,]*betas
length(snp)

#Apply the function to each row (person) to get PRS for each person
all_PRS = apply(snp, 1, compute_PRS)


#Use parallel computing (not done yet)
library(parallel)
no_cores = detectCores()-1
cl - makeCluster(no_cores)

all_PRS <- parSapply(cl=cl, X=snp, FUN=compute_PRS)

stopCluster(cl)