#EDA

## NEXT STEPS:
#
# - Make this script into a function
# - apply to genetic data, which will be much more interesting
# - see if there is way to group thesemeasurements by disease outcome or something
# - investigate use of IDA and LDA
#

setwd("/rdsgpfs/general/user/lc5415/home/hda_tds_ukbiobank")
if (!require(devtools)) install.packages('devtools')
library(devtools)
if (!require(remotes)) install.packages('remotes')
library(remotes)
if (!require(ggbiplot)) install_github("vqv/ggbiplot")
library(ggbiplot)
# Load datasets, when doing actual analysis, replace these by the non-toy datasets
b = readRDS("bio_toy.rds")
c = readRDS("cov_toy.rds")

# Some of b columns are character eventhough they are numeric measures,
# make those into numeric and put back into a df
b = as.data.frame(sapply(b, as.numeric))

# safety-check for all vars being numeric
stopifnot(all(apply(b, 2, is.numeric)))


#remove missing values, this can later be replace by imputation to compare results,
# MAR is probably the approach to take as we are dealing with biochemical measurements
# not humans
b = b[complete.cases(b),]

# compute b's principal components,
# we set center and scale as TRUE, so that all features are scaled (~divided by sd) 
# and centred (~mean removed) before calculating PCAs as it should be.
b.pca = prcomp(b, center = T, scale. = T)

# print summary of pca
summary(b.pca)

# plotting original feature vectors and measurements projected onto top 2 PCAs
ggbiplot(b.pca)
