# Aim of this script is to replicate the work from practical 4 on the biomarker dataset

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(devtools)
if (!require(mixOmics)) devtools::install_github("mixOmicsTeam/mixOmics")
library(sgPLS)
library(pheatmap)

bio <- readRDS("data/preprocessed/bioImputed.rds")
cov <- readRDS("data/preprocessed/covProcessed.rds")

X = bio
y = cov$CVD_status


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

#Sourcing from Barbz code, running calibration to find which number of selected variables gives the lowest misclassification rate. 
source("pls_functions.R")
set.seed(1)
res_splsda = CalibratesPLSDA(X,y, ncomp=1, Nrepeat=10)
PlotCalib(res = res_splsda)


#################################################################
##            Visualising the loadings coefficients            ##
#################################################################

#Creating a plot to visualise the loadings coefficients of all biomarkers 
Loadings = cbind(sPLSDA$loadings$X, rep(NA, 28))
Loadings = as.vector(t(Loadings))
plot(Loadings, col = c('red', NA), xaxt = "n", 
     ylab = "Loadings Coefficients", 
     type = "h", lwd = 3, xlab = "")
axis(1, at = seq(1.5, 28 * 4, by = 4), labels = colnames(PLSDA$X),
     las = 2)
abline(h = 0, lty = 2)




##################################################################
##                      Stability analyses                      ##
##################################################################

#Creating a heatmap of the selection of variables, over 100 iterations each selecting a different training/test set 
set.seed(1)
Stability_results = StabilityPlot(X = X, Y = y, NIter = 100)
pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, filename = "results/PLS_stability.pdf",
         height = 5, width = 10)

