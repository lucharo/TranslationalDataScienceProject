# This file will only include preprocessing of the original datasets and their respective storage

####################################################################
##  Input this file: Original covariates data set (37 columns),   ##
##   Original biomarkers(60 columns), Biomarker annotation file,  ##
##         Covariate dictionary, SNP dataset, SNP info file       ##
####################################################################

#################################################################################
##       Output from this file: Processed covariates, only with entries        ##
##          relevant to CVD[covProcessed.rds], Biomarkers dataset with         ##
##                   well labelled columns and only one entry                  ##
##                  over time (not imputed)[bioProcessed.rds],                 ##
##                     same biomarkers dataset but imputed                     ##
##          through MICE imputation[bioImputed.rds], biomarker dataset         ##
##                   with only complete cases [bioMCAR.rds],                   ##
##   biomarker preprocessed with all columns (30 columns)[bioUnfiltered.rds]   ##
#################################################################################

###########################################################################
###########################################################################
###                                                                     ###
###                         PACKAGE DECLARATION                         ###
###                                                                     ###
###########################################################################
###########################################################################


if (!require(devtools)) install.packages('devtools')
library(devtools)
if (!require(remotes)) install.packages('remotes')
library(remotes)
if (!require(ggbiplot)) install_github("vqv/ggbiplot")
library(ggbiplot)
if (!require(GGally)) install_github("GGally")
library(GGally)
if (!require(tidyverse)) install.packagaes("tidyverse")
library(tidyverse)
if (!require(naniar)) install.packages("naniar")
library(naniar)
if (!require(factoextra)) install.packages("factoextra")
library(factoextra)
if (!require(ggfortify)) install.packages("ggfortify")
library(ggfortify)
if (!require(stats)) install.packages("stats")
library(stats)
if (!require(mice)) install.packages('mice')
library(mice)
if (!require(DMwR)) install.packages('DMwR')
library(DMwR)
if (!require(impute)) BiocManager::install("impute")
library(impute)



library(parallel)
cores = detectCores()

############################################################################
############################################################################
###                                                                      ###
###                             DATA LOADING                             ###
###                                                                      ###
############################################################################
############################################################################



##################################################################
##                        Cluster add-in                        ##
##################################################################
cluster = 1
t0 = Sys.time()
platform = Sys.info()['sysname']
if (cluster == 1 | platform == "Linux"){
  cov.original = readRDS("../FULLDATA/Covariates_full.rds")
  bio.original= readRDS("../FULLDATA/Biomarkers_full.rds")
  bio.dict = readxl::read_xlsx("../Biomarker_annotation.xlsx")
  cov.dict = readxl::read_xlsx("../Covariate_dictionary.xlsx")
  snp.original = readRDS('../FULLDATA/genetic_data_cvd_snps.rds')
  snp_info.original = readxl::read_xlsx("../SNP_info.xlsx")
  
  # id's into columnss, not doing this anymore
  # rownames(bio.original) = bio.original$`mydata$eid`
  # bio.original = bio.original[,-1]
  
  # setting up ID column
  colnames(bio.original)[1] = 'ID'
  cov.original = cbind(ID = rownames(cov.original), cov.original)
  snp.original = cbind(ID = rownames(snp.original), snp.original)
  
  cov.original$ID = as.numeric(levels(cov.original$ID)[cov.original$ID])
  snp.original$ID = as.numeric(levels(snp.original$ID)[snp.original$ID])
  
  save_path = "../FULLDATA/preprocessed/"
} else {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    #################################################################
    ##                      Original datasets                      ##
    #################################################################
    cov.original = readRDS("../data/Covars_toy.rds")
    bio.original= readRDS("../data/Biomarkers_toy.rds")
    bio.dict = readxl::read_xlsx("../Biomarker_annotation.xlsx")
    cov.dict = readxl::read_xlsx("../Covariate_dictionary.xlsx")
    snp.original = readRDS('../data/Genes_toy.rds')
    snp_info.original = readxl::read_xlsx("../SNP_info.xlsx")

    bio.original = cbind(ID = 1:nrow(bio.original), bio.original)
    cov.original = cbind(ID = 1:nrow(cov.original), cov.original)
    snp.original = cbind(ID = 1:nrow(snp.original), snp.original)
    
    save_path = "../data/preprocessed/"
}



##################################################################


##################################################################
##              Changing biomarkers codes by names              ##
##################################################################
#make nicely looking names (programmingly functional)
colnames(bio.dict) = make.names(colnames(bio.dict), unique=TRUE)

#get column numbers of columns with name containing pattern *(.)1(.)*
# use (.) to match the dot as opposed to using . as a wildcard
bio = bio.original[,c(T, !grepl("*(.)1(.)0", colnames(bio.original)[-1]))]

# Match code with biomarker name to change column names of b
# get element 2 to 6 of all string in vector colnames(b)
# the match() function, match the substring from colnames to
# the UK.biobank.field in the biomarkers dictionary, 
# effectively ordering the colnames of b
# Alternative: order UK.bionbank.field entries and match them
#---- bio.dict = bio.dict %>% arrange(UK.Biobank.Field)

colnames(bio)[-1] = bio.dict$Biomarker.name[
  match(substring(colnames(bio)[-1],2,6),bio.dict$UK.Biobank.Field)] 

colnames(bio)[-1] = make.names(colnames(bio)[-1], unique=TRUE)
colnames(bio)[-1] = sub("\\.\\.",".", colnames(bio)[-1])

# safety-check for all vars being numeric
stopifnot(all(apply(bio, 2, is.numeric)))


##################################################################
##               Processing of covariates dataset               ##
##################################################################
## preprocessing c

#remove anything to do wih cancer or external deaths
cov = cov.original[,!grepl("cancer|external", colnames(cov.original))]

#remove codes except for icd10
cov = cov[,!(colnames(cov) %in% c("cvd_final_icd9",
                                  "cvd_final_nc_illness_code",
                                  "cvd_final_opcs4",
                                  "cvd_final_ukb_oper_code",
                                  "other_cause_death"))]

unfactor = function(column, df){
  #' Check if column in integer or string and  turn to the correct type
  #' to avoid other function turning into factor indices
  
  if (class(df[1,column]) == "factor"){
    if (typeof(levels(df[1,column])[df[1,column]]) == "integer"){
      out = as.numeric(levels(df[,column]))[df[,column]]
    } else if (typeof(levels(df[1,column])[df[1,column]])  == "character") {
      out = as.character(levels(df[,column]))[df[,column]]
    } else { #NA
      out = as.character(levels(df[,column]))[df[,column]]
    }
  } else {
    out = df[,column]
  }

  return(out)
}

names.cov = colnames(cov)

cov = as.data.frame(
  lapply(colnames(cov),
       function(column) {
         print(column)
       unfactor(column, cov)}), stringsAsFactors = F
  )

colnames(cov) = names.cov

# replace empty strings by NA
values_to_replace_w_na = c("")
#################################################################
##                     This line is GUILTY                     ##
#################################################################
cov = cov %>% replace_with_na_all(condition = ~.x %in%
                                             values_to_replace_w_na)

cov = cov[!rowMeans(is.na(cov))>0.5,]
saveRDS(cov, file = paste0(save_path,"covProcessed.rds"))


##################################################################
##                   Processing of biomarkers                   ##
##################################################################

saveRDS(bio, file = paste0(save_path,"bioUnfiltered.rds"))
# drop columns with more than 50% missing values i.e. Rheumatoid factor
# and oestradiol
bio = bio[,!colMeans(is.na(bio))>0.5]
# drop rows with more than 50% biomarkers missing
bio = bio[!rowMeans(is.na(bio))>0.5,]
saveRDS(bio, file = paste0(save_path,"bioProcessed.rds"))

bioMCAR = bio[complete.cases(bio),]
saveRDS(bioMCAR, file = paste0(save_path,"bioMCAR.rds"))


#################################################################
##               IMPUTATION MOVED TO OTHER FILES               ##
#################################################################

# impute biomarkers based on biomarkers only
# t0 = Sys.time()
# imp.model = mice(bio, m=5, maxit = 10,
#                  seed = 500, printFlag = F)
# print(Sys.time() - t0) # takes about 1 minute
# 
# # Impute with parallelisation
# t0 = Sys.time()
# imp.model = parlmice(bio,  m =5, seed = NA, printFlag = F,
#                      cl.type = "FORK")
# print(Sys.time() - t0) # takes about 1 minute

# here we assign the imputed data to bio.imp
#bio.imp = complete(imp.model,2)


#################################################################
##                     Proccessing of SNPs                     ##
#################################################################

# removing any snps from the info file that are not included in the data provided by barbz
snps = colnames(snp.original)
snp.info = snp_info.original[snp_info.original$markername %in% snps, ]

saveRDS(snp.info, file = paste0(save_path,"snpInfo.rds"))

snp = cbind(ID = snp.original$ID, 
             as.data.frame(
  apply(snp.original[,-1], 2,
        as.numeric)
  ))

# recoding 00 as NA,  01 as 0, 02 as 1, and 03 as 2
# snp <- snp.original
##################################################################
##           See code above for same thing with apply           ##
##################################################################
# for (i in colnames(snp)){
#   snp[,i] = as.numeric(snp[,i], multiple=TRUE)
# }

snp[snp==0] <- NA
snp[snp==1] <- 0
snp[snp==2] <- 1
snp[snp==3] <- 2

snp$ID = snp.original$ID
# SNP'S DATA A BIT HARDER TO IMPUTE BECAUSE VALUES ARE EITHER 0,1 OR 2, 
# SOME ENTRIES (PARTICUARLY COLUMN 77 HAS A LOT OF THE SAME VALUES, =2 FOR COLUMN 77)
# WHICH MAKES IT EASY TO DO IMPUTATION:
# -- mean(snp[,77]==2, na.rm = T) == 0.995998
# -- sum(colMeans(snp==0, na.rm = T)>0.7) = 0
# -- sum(colMeans(snp==1, na.rm = T)>0.95) = 0
# -- sum(colMeans(snp==2, na.rm = T)>0.95) = 9 # having these columns in results in a standard deviation of
# 0 which does not allow us to scale the SNP values, BUT MAYBE SNPs don't need to be rescaled
# -- sum(colMeans(snp==2, na.rm = T)>0.7) = 60
snp = snp[,!(colMeans(snp==2,na.rm= T)>0.9)]

# NOT REALLY SURE HOW TO IMPUTE SNPs, ASK BARBS

# RMSE.snps = t(sapply(1:10,function(x) kNNImputeOptimization(data.in = snp, seed = x)))
# boxplot(RMSE.snps)
# best.k.med.snps = which.min(apply(RMSE.snps, 2, median))
# best.k.mean.snps = which.min(colMeans(RMSE.snps))

#### FOR NOW, QUICK AND DIRTY SOLUTION
snp.imp = data.frame(impute.knn(as.matrix(snp[,-1]))$data)
snp.imp = cbind(ID = snp$ID, snp.imp)
snp.imp = as.data.frame(snp.imp, stringsAsFactors = FALSE)
saveRDS(snp.imp, paste0(save_path,"snpImputed.rds"))

saveRDS(snp, paste0(save_path,"snpProcessed.rds"))

###Can do imputation of snps
# First check that no people have all missing values and 
# no snps have all missing values 
max(rowSums(is.na(snp)))
#max number of missing values for one person is 15/177
max(colSums(is.na(snp)))
#max number of missing values for one snp is 226/2000

# list <-lapply(1:10,
#               function(col) ggplot2::qplot(snp.imp[[col]],
#                                            geom = "histogram",
#                                            binwidth = 1))
# 
# cowplot::plot_grid(plotlist = list)
print(paste0("Common rows between bio and cov: ",sum(bio$ID %in% cov$ID)))
print(Sys.time()-t0)
      