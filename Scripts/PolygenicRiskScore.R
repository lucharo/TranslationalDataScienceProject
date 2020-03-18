# Aim of this script is to calculate a polygenic risk score for each person

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

snp = readRDS(paste0(data_folder,"snpImputed.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
snp.info = readRDS(paste0(data_folder,"snpInfo.rds"))


#Correcting snp IDs using Barbz bridge file
bf = readRDS(paste0(data_folder,"bridge_file.rds"))
bf = as.data.frame(bf)
snp.new = merge(snp, bf, by.x='ID', by.y='previous')
snp.new = snp.new[,-1]
colnames(snp.new)[163] = 'ID'


##################################################################
##                        Computing PRS                         ##
##################################################################

# Added this for cluster as snp file is smaller in cluster than ni toy? 162 snps vs 177?
snps = colnames(snp.new)
snp.info = snp.info[snp.info$markername %in% snps, ]

#Check that the snps are in the same order in snp and snp.info ([-1] to remove the ID column)
stopifnot(all(colnames(snp.new)[-163] == snp.info$markername))

#Extracting the beta values (which are in the correct order based on above check)
betas <- snp.info$beta

#Element-wise multiplication of each row of snp by the betas
# (margin=2 specifies it is rows of the matrix; each row is a person)
# this multiples each no. of snp copies by the beta coefficient for 
# that snp 
PRS_matrix = sweep(snp.new[,-163], MARGIN=2, betas, `*`)


#Now calculate the weighted sum for each person (have included na.rm=TRUE to make it work for now, can remove when we've done imputation)
PRS_sums = rowSums(PRS_matrix, na.rm=TRUE)

PRS = cbind(snp.new$ID, PRS_sums)
PRS = as.data.frame(PRS)
colnames(PRS) = c("ID", "PRS")
saveRDS(PRS, paste0(save_data,"PolygenicRiskScore.rds"))


#################################################################
##                Comparison between CVD status                ##
#################################################################

## Merge cov and PRS
cov.prs <- merge(cov, PRS, by="ID")

#Boxplot of PRS by CVD status
prs_boxplot <- ggplot(cov.prs, aes(x=CVD_status, y=PRS)) +
  geom_boxplot()
saveRDS(prs_boxplot, paste0(save_plots,"PRS_boxplot.rds"))

#t-test - no sig difference in mean PRS between groups...
t_test = t.test(PRS ~ CVD_status, data=cov.prs)
saveRDS(t_test, paste0(save_plots,"PRS_ttest.rds"))




