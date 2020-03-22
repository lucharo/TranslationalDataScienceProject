# Aim of this script is to run the sPLS-DA and sgPLS-DA models with calibrated parameters

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())

suppressPackageStartupMessages(library(devtools))
#if (!require(mixOmics)) devtools::install_github("mixOmicsTeam/mixOmics")
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(sgPLS))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))


cluster = 0

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

bio <- readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov <- readRDS(paste0(data_folder,"covProcessed.rds"))

cvd <- cov %>% select(ID, CVD_status)
bio.cov <- merge(bio, cvd, by='ID')


##################################################################
##                   Fitting sPLS-DA models                     ##
##################################################################

#Select all biomarkers from bio.cov for X
X = bio.cov[, 2:29]
y = bio.cov$CVD_status


#Sparse plsda model; keepX is number of parameters to keep (9 had lowest misclassification rate in calibration). The final line returns the variables selected.
sPLSDA <- splsda(X, y, ncomp=1, mode='regression', keepX=9)
sPLSDA$loadings
sPLSDA$explained_variance
sPLSDA$loadings$X[sPLSDA$loadings$X != 0, ]

##Since the loading for y(1) is positive then the interpretation of the x loadings is as normal (positive loadings are higher in cases).  



