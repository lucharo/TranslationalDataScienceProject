# Aim of this script is to replicate the work from practical 4 on the biomarker dataset

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())

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

bio <- readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov <- readRDS(paste0(data_folder,"covProcessed.rds"))

cvd <- cov %>% select(ID, CVD_status)
bio.cov <- merge(bio, cvd, by='ID')


#Create training and test sets 
smp_size <- floor(0.8*nrow(bio.cov))

#set the seed to make the partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(bio.cov)), size = smp_size)

train <- bio.cov[train_ind, ]
test <- bio.cov[-train_ind, ]

#Select all biomarkers from bio.cov for X
X_train = train[, 2:29]
X_test = test[, 2:29]
y_train = train$CVD_status
y_test = test$CVD_status


#################################################################
##                      Model calibration                      ##
#################################################################

#Sourcing from Barbz code, this function performs 5-fold cross-validation to find which number of selected variables gives the lowest misclassification rate (of CVD status). 
source("pls_functions.R")
set.seed(1)
res_splsda = CalibratesPLSDA(X_train, y_train, ncomp=1, Nrepeat=100)

pdf(paste0(save_plots,"sPLSDA_calibration.pdf"))
splsda_calibration <- PlotCalib(res = res_splsda)
dev.off()
saveRDS(splsda_calibration, paste0(save_plots,"sPLSDA_calibration.rds"))



##################################################################
##                       Sparse group pls                       ##
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
  

#Making the gPLSDA model with fran's groupings. Group membership is set with ind.block.x. 
#The number of groups to be selected by the model is defined by keepX. 
gPLSDA <- gPLSda(X_fran, y, ncomp = 1, 
                 ind.block.x = X_cuts_fran, keepX = 1)

gPLSDA$loadings$X
#can see that the endocrine group is selected.

#Re-run with paper's groups:
gPLSDA <- gPLSda(X_paper, y, ncomp = 1, 
                 ind.block.x = X_cuts_paper, keepX = 1)

gPLSDA$loadings$X
#This time the liver group is selected. If you run fran's model now with keepX=2, then liver is next selected. 

#Sparse group plsda introduces another layer of sparsity - selecting variables within the group. 
#This is done with alpha.x (between 0 and 1; the closer it is to 1, the more coefficients are shrunk)

sgPLSDA <- sgPLSda(X_fran, y, ncomp = 1,
                  ind.block.x = X_cuts_fran, keepX = 1, alpha.x = 0.1)
sgPLSDA$loadings$X



#################################################################
##               Model calibration for group PLS               ##
#################################################################

#The two parameters of the sgPLS-DA model can be calibrated using the function CalibratesgPLSDA()
set.seed(1)
res_sgplsda = CalibratesgPLSDA(dataX = X_fran, dataY = y, ncomp = 1,
                               Nrepeat = 50, Xgroups = X_cuts_fran)
pdf(paste0(save_plots,"sgPLSDA_calibration.rds"))
sgplsda_calibration <- PlotCalib(res = res_sgplsda, type = "sgPLSDA")
saveRDS(sgplsda_calibration, paste0(save_plots,"sgPLSDA_calibration.rds"))


