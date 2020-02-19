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
  PRS = snp[i,]*betas[i]
}


compute_PRS <- function(){
  PRS = c()
  for (i in length(snp)){
    PRS[i] = snp[i,]*betas[i]
  }
}


# Using apply function to compute PRS for all subjects 
all_PRS <- apply(snp, MARGIN=1, FUN=compute_PRS)
