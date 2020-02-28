# Aim of this script is to replicate the work from practical 4 on the biomarker dataset

##################################################################
##                 Prepare libraries and data                   ##
##################################################################


## PROBLEM WITH PLS NOT WORKING HAS TO DO WITH CASE-CONTROLSIMBALANCE,
# OUT OF 2000 LIKE 73 CASES OR SO

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(devtools)
if (!require(mixOmics)) devtools::install_github("mixOmicsTeam/mixOmics")
library(sgPLS)
library(pheatmap)

bio <- readRDS("data/preprocessed/bioImputed.rds")
cov <- readRDS("data/preprocessed/covProcessed.rds")

bio.cov <- merge(bio, cov['CVD_status'], by='row.names')
X = as.data.frame(bio)
y = bio.cov$CVD_status

##################################################################
##                    Running PLS-DA model                    ##
##################################################################

#Non-penalised plsda (i.e. no feature selection)
PLSDA <- plsda(X, y, ncomp=1, mode="regression")
PLSDA$loadings
PLSDA$explained_variance

#Sparse plsda model; keepX is number of parameters to keep. The final line returns the variables selected.
sPLSDA <- splsda(X, y, ncomp=1, mode='regression', keepX=5)
sPLSDA$loadings
sPLSDA$explained_variance
sPLSDA$loadings$X[sPLSDA$loadings$X != 0, ]


#################################################################
##                      Model calibration                      ##
#################################################################

#Sourcing from Barbz code, this function performs 5-fold cross-validation to find which number of selected variables gives the lowest misclassification rate (of CVD status). 
source("pls_functions.R")
set.seed(1)
res_splsda = CalibratesPLSDA(X,y, ncomp=1, Nrepeat=10)
PlotCalib(res = res_splsda)

#Constant misclassification rate, which is the same as the proportion of CVD cases in the dataset:
sum(cov$CVD_status == 2)/length(cov$CVD_status)
#i.e. the PLS regression does not predict any CVD cases




##################################################################
##                      Stability analyses                      ##
##################################################################

#Creating a heatmap of the selection of variables, over 100 iterations each selecting a different training/test set 
set.seed(1)
Stability_results = StabilityPlot(X = X, Y = y, NIter = 100)
pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, filename = "results/PLS_stability.pdf",
         height = 5, width = 10)



##################################################################
##                       Trying group pls                       ##
##################################################################

#List of biomarkers in order of groups (based on Fran's grouping): first 8 are liver, the next 10 are metabolic, next 2 immune, next 5 endocrine, and final 3 kidney. X_cuts defines these cuts. 
groups_fran = c('Alanine.aminotransferase','Alkaline.phosphatase','Aspartate.aminotransferase','Direct.bilirubin','Gamma.glutamyltransferase','Total.bilirubin','Total.protein','Albumin','Apolipoprotein.A','Apolipoprotein.B','Cholesterol','Glucose','Glycated.haemoglobin.HbA1c.','HDL.cholesterol','LDL.direct','Lipoprotein.A','Triglycerides','Urate','C.reactive.protein','IGF.1','Calcium','Phosphate','SHBG','Testosterone','Vitamin.D','Creatinine','Cystatin.C','Urea')

X_fran = X[, groups_fran]
X_cuts_fran = c(8, 18, 20, 25)

#List of biomarkers in order of groups (based on paper's grouping): first 3 are liver, the next 4 are metabolic, next 2 immune, next 1 endocrine, and final 1 kidney (this grouping is the same but omitting some). X_cuts again defines the cuts.
groups_paper = c('Alanine.aminotransferase','Aspartate.aminotransferase','Gamma.glutamyltransferase','Cholesterol','Glycated.haemoglobin.HbA1c.','HDL.cholesterol','Triglycerides','C.reactive.protein','IGF.1','Testosterone','Creatinine')

X_paper = X[, groups_paper]
X_cuts_paper = c(3, 7, 9, 10)
  


#Making the gPLSDA model with fran's groupings. Group membership is set with ind.block.x. The number of groups to be selected by the model is defined by keepX. 
gPLSDA <- gPLSda(X_fran, y, ncomp = 1, 
                 ind.block.x = X_cuts_fran, keepX = 1)

gPLSDA$loadings$X
#can see that the endocrine group is selected.

#Re-run with paper's groups:
gPLSDA <- gPLSda(X_paper, y, ncomp = 1, 
                 ind.block.x = X_cuts_paper, keepX = 1)

gPLSDA$loadings$X
#This time the liver group is selected. If you run fran's model now with keepX=2, then liver is next selected. 

#Sparse group plsda introduces another layer of sparsity - selecting variables within the group. This is done with alpha.x (between 0 and 1; the closer it is to 1, the more coefficients are shrunk)

sgPLSDA <- sgPLSda(X_fran, y, ncomp = 1,
                  ind.block.x = X_cuts_fran, keepX = 1, alpha.x = 0.1)
sgPLSDA$loadings$X


#################################################################
##               Model calibration for group PLS               ##
#################################################################

#The two parameters of the sgPLS-DA model can be calibrated using the function CalibratesgPLSDA()
set.seed(1)
res_sgplsda = CalibratesgPLSDA(dataX = X_fran, dataY = y, ncomp = 1,
                               Nrepeat = 5, Xgroups = X_cuts_fran)
PlotCalib(res = res_sgplsda, type = "sgPLSDA")



#################################################################
##            Visualising the loadings coefficients            ##
#################################################################

#################################################################
##            Visualising the loadings coefficients            ##
#################################################################

#This plot visualises the loadings coefficients of all biomarkers, using the sPLDSA model created earlier 

Loadings = cbind(sPLSDA$loadings$X, rep(NA, 28))
Loadings = as.vector(t(Loadings))
plot(Loadings, col = c('red', NA), xaxt = "n", 
     ylab = "Loadings Coefficients", 
     type = "h", lwd = 3, xlab = "")
axis(1, at = seq(1.5, 28 * 4, by = 4), labels = colnames(PLSDA$X),
     las = 2)
abline(h = 0, lty = 2)


#This plot visualises the loadings coefficients obtained from both sPLSDA and sgPLSDA models

Loadings = cbind(sPLSDA$loadings$X, sgPLSDA$loadings$X,
                 rep(NA, 28), rep(NA, 28))
Loadings = as.vector(t(Loadings))
Loadings = Loadings[-c(length(Loadings) - 1, length(Loadings))]
par(mar = c(10, 5, 3, 3))
plot(Loadings, col = c('red', 'blue', NA, NA),
     xaxt = "n", ylab = "Loadings Coefficients", type = "h",
     lwd = 3, xlab = "")
axis(1, at = seq(1.5, 28 * 4, by = 4), labels = colnames(PLSDA$X), las = 2)
axis(1, at = c(0, X_cuts_fran, 28) * 4, line = 6, labels = NA)
#axis(1, at = c(X_cuts_fran) * 4, 
 #    labels = c("Liver","Metabolic", "Immune","Endocrine","Kidney"),
  #   line = 6, tick = FALSE)
abline(v = c(0, Xgroups, 28) * 4, lty = 3, col = "black")
abline(h = 0, lty = 2)
legend("topleft", legend = c("sPLS-DA", "sgPLS-DA"),
       lty = 1, lwd = 3, col = c('red', 'blue'), cex = 0.5)
