args = commandArgs(trailingOnly = TRUE)
seed = as.numeric(args[1])

#Aim of this script is to calibrate the sparse group PLS-DA model (parameters = number of groups to include and sparsity parameter)              

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

#rm(list=ls())

suppressPackageStartupMessages(library(devtools))
#if (!require(mixOmics)) devtools::install_github("mixOmicsTeam/mixOmics")
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(sgPLS))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))


cluster = 1

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

ifelse(!dir.exists(file.path(save_plots, "PLS/")),
       dir.create(file.path(save_plots, "PLS/")), FALSE)
save_plots = paste0(save_plots,"PLS/")

bio <- readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov <- readRDS(paste0(data_folder,"covProcessed.rds"))

cvd <- cov %>% select(ID, CVD_status)
bio.cov <- merge(bio, cvd, by='ID')

#Select all biomarkers from bio.cov for X
X = bio.cov[, 2:29]
y = bio.cov$CVD_status


##################################################################
##                Setting up groups for sgPLS                   ##
##################################################################

#List of biomarkers in order of groups (based on Fran's grouping): first 8 are liver, the next 10 are metabolic, next 2 immune, next 5 endocrine, and final 3 kidney. X_cuts defines these cuts. 
groups_fran = c('Alanine.aminotransferase','Alkaline.phosphatase','Aspartate.aminotransferase',
                'Direct.bilirubin','Gamma.glutamyltransferase','Total.bilirubin',
                'Total.protein','Albumin','Apolipoprotein.A','Apolipoprotein.B','Cholesterol',
                'Glucose','Glycated.haemoglobin.HbA1c.','HDL.cholesterol','LDL.direct',
                'Lipoprotein.A','Triglycerides','Urate','C.reactive.protein','IGF.1','Calcium',
                'Phosphate','SHBG','Testosterone','Vitamin.D','Creatinine','Cystatin.C','Urea')

X_fran = X[, groups_fran]
X_cuts_fran = c(8, 18, 20, 25)

#List of biomarkers in order of groups (based on paper's grouping): first 3 are liver, the next 4 are metabolic, next 2 immune, next 1 endocrine, and final 1 kidney (this grouping is the same but omitting some). X_cuts again defines the cuts.
groups_paper = c('Alanine.aminotransferase','Aspartate.aminotransferase',
                 'Gamma.glutamyltransferase','Cholesterol','Glycated.haemoglobin.HbA1c.',
                 'HDL.cholesterol','Triglycerides','C.reactive.protein','IGF.1','Testosterone',
                 'Creatinine')

X_paper = X[, groups_paper]
X_cuts_paper = c(3, 7, 9, 10)

  
print("Data loaded")
#################################################################
##                 Model calibration for sgPLS                 ##
#################################################################

#Sourcing from Barbz code
#The two parameters of the sgPLS-DA model can be calibrated using CalibratesgPLSDA()
#This function performs 5-fold cross-validation to find which paramater combination gives the lowest misclassification rate (of CVD status). 

source("pls_functions.R")
seed = as.numeric(args[1])
set.seed(seed)
print("PLS started")
res_sgplsda = CalibratesgPLSDA(dataX = X_fran, dataY = y, ncomp = 1,
                               Nrepeat = 1, Xgroups = X_cuts_fran)
print("PLS finished")

ifelse(!dir.exists(file.path(save_plots, "ArrayJob/")),
       dir.create(file.path(save_plots, "ArrayJob/")), FALSE)
save_plots = paste0(save_plots,"ArrayJob/")

saveRDS(res_sgplsda, paste0(save_plots,"sgCalibration_",as.character(seed),".rds"))

#pdf(paste0(save_plots,"sgPLSDA_calibration.rds"))
#sgplsda_calibration <- PlotCalib(res = res_sgplsda, type = "sgPLSDA")
#saveRDS(sgplsda_calibration, paste0(save_plots,"sgPLSDA_calibration.rds"))


