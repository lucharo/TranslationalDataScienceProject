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

#Check that the snps are in the same order in snp and snp.info
stopifnot(all(colnames(snp) == snp.info$markername))

#Extracting the beta values (which are in the correct order based on above check)
betas <- snp.info$beta

#Element-wise multiplication of each row of snp by the betas (margin=2 specifies it is rows of the matrix; each row is a person) - this multiples each no. of snp copies by the beta coefficient for that snp 
PRS_int = sweep(snp, MARGIN=2, betas, `*`)


#Now calculate the weighted sum for each person (have included na.rm=TRUE to make it work for now, can remove when we've done imputation)
PRS = rowSums(PRS_int, na.rm=TRUE)


#################################################################
##                Comparison between CVD status                ##
#################################################################

cov <- readRDS("data/preprocessed/covProcessed.rds")
cov.prs <- cbind(cov, PRS)

#Boxplot of PRS by CVD status
ggplot(cov.prs, aes(x=CVD_status, y=PRS)) + geom_boxplot()

#t-test - no sig difference in mean PRS between groups...
t.test(PRS ~ CVD_status, data=cov.prs)
