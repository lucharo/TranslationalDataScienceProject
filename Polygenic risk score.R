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

# Function for computing PRS per person
betas <- snp.info$beta

# establishing for loop to use inside function 
for (i in length(snp)){
  for (j in length(snp)){
    PRS = snp[i,j]*betas
  }
}

PRSt = snp[1,]*betas
PRSx = rowSums(PRSt)

compute_PRS <- function(){
  sb = c()
  for (i in length(snp)){
    sb[i,] = snp[i,]*betas
  }
  PRS = rowSums(sb)
}

##this for loop calculates beta*copy number for each snp  

snps.betas <- function() {
  for (i in length(snp)){
    PRSt = snp[i,]*betas
  }
}

snps.betas2 <- function(row = length(snp)){
  PRSt = snp['row',]*betas
}

PRSt = apply(snp, 1, snps.betas)

compute_PRS <- function(){
  sums_rows = apply(snp, 1, sums)
  sb = snp[i,]*beta
  PRS = rowSums(sb)
}


# Using apply function to compute PRS for all subjects 
all_PRS <- apply(snp, MARGIN=1, FUN=compute_PRS)
