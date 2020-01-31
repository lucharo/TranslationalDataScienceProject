#EDA

## NEXT STEPS:
#
# - Make this script into a function
# - apply to genetic data, which will be much more interesting
# - see if there is way to group thesemeasurements by disease outcome or something
# - investigate use of IDA and LDA
#

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!require(devtools)) install.packages('devtools')
library(devtools)
if (!require(remotes)) install.packages('remotes')
library(remotes)
if (!require(ggbiplot)) install_github("vqv/ggbiplot")
library(ggbiplot)
if (!require(GGally)) install_github("GGally")
library(GGally)
if (!require(tidyverse)) install("tidyverse")
library(tidyverse)
if (!require(naniar)) install("naniar")
library(naniar)


## LOAD DATASETS ####
# Load datasets, when doing actual analysis, replace these by the non-toy datasets
cov.original = readRDS("data/Covars_toy.rds")
bio.original= readRDS("data/Biomarkers_toy_example.rds")
bio.dict = readxl::read_xlsx("Biomarker_annotation.xlsx")

cov.dict = readxl::read_xlsx("Covariate_dictionary.xlsx")

## Initial pre-processing ####
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

# check number of missing values per column
colSums(is.na(bio))

# safety-check for all vars being numeric
stopifnot(all(apply(bio, 2, is.numeric)))

## preprocessing c
# replace empty strings by NA
values_to_replace_w_na = c("")
cov = cov.original %>% replace_with_na_all(condition = ~.x %in% values_to_replace_w_na)
#remove anything to do wih cancer or external deaths
cov = cov[,!grepl("cancer|external", colnames(cov))]

## Second pre-processing ####

# merging b with c
bio.joint = merge(bio,cov,by="row.names",all.x=TRUE)[,c(colnames(bio),"vit_status")]

#remove missing values, this can later be replace by imputation to compare results,
# MAR is probably the approach to take as we are dealing with biochemical measurements
# not humans
bio.joint = bio.joint[complete.cases(bio.joint),]

# compute b's principal components,
# we set center and scale as TRUE, so that all features are scaled (~divided by sd) 
# and centred (~mean removed) before calculating PCAs as it should be.
b.pca = prcomp(bio, center = T, scale. = T)

# print summary of pca
summary(b.pca)

# plotting original feature vectors and measurements projected onto top 2 PCAs
ggbiplot(b.pca)


#################### Catriona having fun

plot(c$age_cancer,c$age_recruitment.0.0)

ggpairs(bio)

