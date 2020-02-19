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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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

library(parallel)
cores = detectCores()

############################################################################
############################################################################
###                                                                      ###
###                             DATA LOADING                             ###
###                                                                      ###
############################################################################
############################################################################

#################################################################
##                      Original datasets                      ##
#################################################################
cov.original = readRDS("data/Covars_toy.rds")
bio.original= readRDS("data/Biomarkers_toy.rds")
bio.dict = readxl::read_xlsx("Biomarker_annotation.xlsx")

##################################################################
##                        Cluster add-in                        ##
##################################################################
cluster = 0

if (cluster == 1){
  rownames(bio.original) = bio.original$`mydata$eid`
  bio.original = bio.original[,-1]
}

##################################################################
cov.dict = readxl::read_xlsx("Covariate_dictionary.xlsx")
snp.original = readRDS('data/Genes_toy.rds')
snp_info.original = readxl::read_xlsx("SNP_info.xlsx")


##################################################################
##              Changing biomarkers codes by names              ##
##################################################################
#make nicely looking names (programmingly functional)
colnames(bio.dict) = make.names(colnames(bio.dict), unique=TRUE)

#get column numbers of columns with name containing pattern *(.)1(.)*
# use (.) to match the dot as opposed to using . as a wildcard
bio = bio.original[,!grepl("*(.)1(.)0", colnames(bio.original))]

# Match code with biomarker name to change column names of b
# get element 2 to 6 of all string in vector colnames(b)
# the match() function, match the substring from colnames to
# the UK.biobank.field in the biomarkers dictionary, 
# effectively ordering the colnames of b
# Alternative: order UK.bionbank.field entries and match them
#---- bio.dict = bio.dict %>% arrange(UK.Biobank.Field)

colnames(bio) = bio.dict$Biomarker.name[
  match(substring(colnames(bio),2,6),bio.dict$UK.Biobank.Field)] 

colnames(bio) = make.names(colnames(bio), unique=TRUE)
colnames(bio) = sub("\\.\\.",".", colnames(bio))

# safety-check for all vars being numeric
stopifnot(all(apply(bio, 2, is.numeric)))


##################################################################
##               Processing of covariates dataset               ##
##################################################################
## preprocessing c
# replace empty strings by NA
values_to_replace_w_na = c("")
cov = cov.original %>% replace_with_na_all(condition = ~.x %in%
                                             values_to_replace_w_na)
#remove anything to do wih cancer or external deaths
cov = cov[,!grepl("cancer|external", colnames(cov))]
#remove codes except for icd10
cov = cov[,!(colnames(cov) %in% c("cvd_final_icd9",
                                  "cvd_final_nc_illness_code",
                                  "cvd_final_opcs4",
                                  "cvd_final_ukb_oper_code",
                                  "other_cause_death"))]

# change numerical binary outcome variables to categorical
str(cov)
cov$CVD_status = as.factor(cov$CVD_status)
cov$vit_status = as.factor(cov$vit_status)
cov$dc_cvd_st = as.factor(cov$dc_cvd_st)
cov$cvd_death = as.factor(cov$cvd_death)

saveRDS(cov, file = "data/preprocessed/covProcessed.rds")


##################################################################
##                   Processing of biomarkers                   ##
##################################################################

saveRDS(bio, file = "data/preprocessed/bioUnfiltered.rds")
# drop columns with more than 50% missing values i.e. Rheumatoid factor
# and oestradiol
bio = bio[,!colMeans(is.na(bio))>0.5]
saveRDS(bio, file = "data/preprocessed/bioProcessed.rds")

bioMCAR = bio[complete.cases(bio),]
saveRDS(bioMCAR, file = "data/preprocessed/bioMCAR.rds")

# impute biomarkers based on biomarkers only
t0 = Sys.time()
imp.model = mice(bio, m=5, maxit = 10,
                 seed = 500, printFlag = F)
print(Sys.time() - t0) # takes about 1 minute

# Impute with parallelisation
t0 = Sys.time()
imp.model = parlmice(bio, n.core = cores, m =5,
                     #m=5, seed = NA, printFlag = F,
                     cl.type = "FORK")
print(Sys.time() - t0) # takes about 1 minute

# here we assign the imputed data to bio.imp
bio.imp = complete(imp.model,2)
saveRDS(bio.imp, file = "data/preprocessed/bioImputed.rds")



#################################################################
##                     Proccessing of SNPs                     ##
#################################################################

# removing any snps from the info file that are not included in the data provided by barbz
snps = colnames(snp.original)
snp.info = snp_info.original[snp_info.original$markername %in% snps, ]

saveRDS(snp.info, file = "data/preprocessed/snpInfo.rds")


# recoding 00 as NA,  01 as 0, 02 as 1, and 03 as 2
snp <- snp.original

for (i in colnames(snp)){
  snp[,i] = as.numeric(snp[,i], multiple=TRUE)
}

snp[snp==0] <- NA
snp[snp==1] <- 0
snp[snp==2] <- 1
snp[snp==3] <- 2

saveRDS(snp, "data/preprocessed/snpProcessed.rds")

