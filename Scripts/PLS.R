# Aim of this script is to replicate the work from practical 4 on the biomarker dataset

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

## PROBLEM WITH PLS NOT WORKING HAS TO DO WITH CASE-CONTROLSIMBALANCE,
# OUT OF 2000 LIKE 73 CASES OR SO

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


##################################################################
##               Basic PLS-DA and sPLS-DA models                ##
##################################################################

#Select all biomarkers from bio.cov for X
X = bio.cov[, 2:29]
y = bio.cov$CVD_status

#Non-penalised plsda (i.e. no feature selection)
PLSDA <- plsda(X, y, ncomp=1, mode="regression")
#PLSDA$loadings
#PLSDA$explained_variance

#Sparse plsda model; keepX is number of parameters to keep. The final line returns the variables selected.
sPLSDA <- splsda(X, y, ncomp=1, mode='regression', keepX=5)
#sPLSDA$loadings
#sPLSDA$explained_variance
#sPLSDA$loadings$X[sPLSDA$loadings$X != 0, ]



#################################################################
##                      Model calibration                      ##
#################################################################

#Sourcing from Barbz code, this function performs 5-fold cross-validation to find which number of selected variables gives the lowest misclassification rate (of CVD status). 
source("pls_functions.R")
set.seed(1)
res_splsda = CalibratesPLSDA(X, y, ncomp=1, Nrepeat=100)

#Make parallel
#no_cores=detectCores()-1
#cl <- makeCluster(no_cores) 
#res_splsda = parSapply(cl = cl, X = X,
#                       FUN = CalibratesPLSDA(X, y, ncomp=1, Nrepeat=100)) 

pdf(paste0(save_plots,"sPLSDA_calibration.pdf"))
splsda_calibration <- PlotCalib(res = res_splsda)
dev.off()
saveRDS(splsda_calibration, paste0(save_plots,"sPLSDA_calibration.rds"))


#Constant misclassification rate, which is the same as the proportion of CVD cases in the dataset:
#sum(cov$CVD_status == 2)/length(cov$CVD_status)
#i.e. the PLS regression does not predict any CVD cases



##################################################################
##                      Stability analyses                      ##
##################################################################

#Creating a heatmap of the selection of variables, over 100 iterations each selecting a different training/test set 
set.seed(1)
Stability_results = StabilityPlot(X = X, Y = y, NIter = 100)
pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, 
         filename = paste0(save_plots,"PLS_stability.pdf"),
         height = 5, width = 10)


##################################################################
##                       Sparse group pls                       ##
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
                               Nrepeat = 50, Xgroups = X_cuts_fran)
pdf(paste0(save_plots,"sgPLSDA_calibration.rds"))
sgplsda_calibration <- PlotCalib(res = res_sgplsda, type = "sgPLSDA")
saveRDS(sgplsda_calibration, paste0(save_plots,"sgPLSDA_calibration.rds"))


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

results = data.frame(rbind(
  cbind(Biomarker = colnames(X),
        Model = 'sPLSDA',
        Loadings = sPLSDA$loadings$X),
  cbind(Biomarker = colnames(X),
        Model = 'sgPLSDA',
        Loadings = sgPLSDA$loadings$X)
))

results = results %>%
  mutate(belong_to = ifelse(Biomarker %in% groups_fran[1:8], "Liver",
                            ifelse(Biomarker %in% groups_fran[9:18], "Metabolic",
                                   ifelse(Biomarker %in% groups_fran[19:20], "Immune",
                                          ifelse(Biomarker %in% groups_fran[21:25], "Endocrine",
                                                 "Kidney")))))

colnames(results)[3] = 'Loadings'
results$minLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) min(0, x)))

results$maxLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) max(0, x)))

PLSDA_loadings = results %>% ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
                              ymax = maxLoad, color = Model))+
  geom_linerange(stat = "identity", position = position_dodge(0.9))+
  geom_point(aes(y = 0), position = position_dodge(0.9)) +
  ylab("Loading coefficients") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(cols = vars(belong_to), scales = "free", space = "free_x")

ggsave(paste0(save_plots,"PLSDA_loadings.pdf"), plot=PLSDA_loadings)
saveRDS(PLSDA_loadings, paste0(save_plots,"PLSDA_loadings.rds"))



#################################################################
##                     Stratified analyses                     ##
#################################################################

#Aim: apply PLS models on subsets of the data with all controls and cases of one particular subtype only. From full data, selected CVD subtypes (based on ICD10 code for death) with more than 100 cases: G454, G459, I200, I209, I210, I211, I214, I219, I249, I251, I259, I635, I639 and I64. 

cvd_icd <- cov %>% select(ID, cvd_final_icd10)
bio.icd <- merge(bio, cvd_icd, by='ID')
y = bio.cov$CVD_status

#Create the stratified datasets 
for (subtype in c("G454", "G459", "I200", "I209", "I210", "I211", "I214",
                  "I219", "I249", "I251", "I259", "I635", "I639", "I64")) {
  Xtemp = X[X$cvd_final_icd10 %in% subtype, ]
  Ytemp = y[X$cvd_final_icd10 %in% subtype]
  assign(paste0("X_", subtype), Xtemp)
  assign(paste0("Y_", subtype), Ytemp)
}

#Computing the misclassification rate by subtype of CVD, for the sPLS-DA and sgPLS-DA models
y_pred <- predict(sPLSDA, newdata = X)
fitted = y_pred$class$max.dist
table(fitted)

MSEP_sPLSDA = NULL
for (subtype in c("","G454","G459","I200","I209","I210","I211","I214",
                  "I219","I249","I251","I259","I635","I639","I64")) {
  idx = which(cov$cvd_final_icd10 == subtype)
  MSEP_sPLSDA[[subtype]] = 1 - sum(diag(table(y[idx],
                               y_pred$class$max.dist[idx])))/length(idx)
}

y_pred_g = predict(sgPLSDA, newdata = X)
fitted_g = y_pred_g$class$max.dist
MSEP_sgPLSDA = NULL
for (subtype in c("","G454","G459","I200","I209","I210","I211","I214",
                  "I219","I249","I251","I259","I635","I639","I64")) {
  idx = which(cov$cvd_final_icd10 == subtype)
  MSEP_sgPLSDA[[subtype]] = 1 - sum(diag(table(y[idx],
                               y_pred_g$class$max.dist[idx])))/length(idx)
}


#Creating a plot of these misclassification rates by CVD subtype

mis_rate = data.frame(rbind(
  cbind(Prediction=fitted, model="sPLSDA", 
        subtype = bio.icd$cvd_final_icd10, truth=bio.cov$CVD_status),
  cbind(Prediction=fitted_g, model="sgPLSDA", 
        subtype = bio.icd$cvd_final_icd10, truth=bio.cov$CVD_status)))

levels(mis_rate$subtype) = c("G454","G459","I200","I209","I210","I211","I214",
                             "I219","I251","I259","I638","I639","I64","I679",
                             "Control")
mis_rate$subtype = replace(mis_rate$subtype, which(is.na(mis_rate$subtype)),
                           "Control")

mis_rate$IncorrectClass = !(mis_rate$comp1 == mis_rate$truth)

mis_rate_plot <- mis_rate %>% filter(subtype %in% c("Control","G454","G459",
                              "I200","I209","I210","I211","I214","I219","I249",
                              "I251","I259","I635","I639","I64")) %>% 
  group_by(subtype, model) %>% 
  summarise(rate = mean(IncorrectClass)) %>% 
  arrange(subtype)

plsda_stratified <- mis_rate_plot %>% ggplot(aes(x = subtype, ymin = 0, ymax = rate, color = model)) + 
  geom_linerange(stat = "identity", position = position_dodge(0.9)) + 
  scale_color_brewer(palette = "Set1") + 
  ylab("Misclassification Rate")

ggsave(paste0(save_plots,"PLSDA_stratified.pdf"), plot=plsda_stratified)
saveRDS(plsda_stratified, paste0(save_plots,"PLSDA_stratified.rds"))


##################################################################
##                      Fitting PLS models                      ##
##################################################################

#Next step is to fit the sPLS-DA, sgPLS-DA and stratified PLS-DA models using the best parameters found earlier by calibration (which we can't do yet until we've run this on the full dataset?)
