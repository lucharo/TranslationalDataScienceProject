# Aim of this script is to run the sPLS-DA and sgPLS-DA models with calibrated parameters

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())

suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(sgPLS))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))


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
##                     Basic PLS-DA model                       ##
##################################################################

#Select all biomarkers from bio.cov for X
X = bio.cov[, 2:29]
y = bio.cov$CVD_status

#Non-penalised plsda (i.e. no feature selection)
PLSDA <- plsda(X, y, ncomp=1, mode="regression")
PLSDA$loadings
PLSDA$explained_variance

##Since the loading for y(1) is positive then the interpretation of the x loadings is as normal 
#(positive loadings are higher in cases).  

#This plot visualises the loadings coefficients obtained from this PLSDA model
results = data.frame(cbind(Biomarker = colnames(X), Loadings = PLSDA$loadings$X))

colnames(results)[2] = 'Loadings'
results$minLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) min(0, x)))
results$maxLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) max(0, x)))

PLSDA_loadings = results %>% ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
                                        ymax = maxLoad))+
  geom_linerange(stat = "identity", position = position_dodge(0.9))+
  geom_point(aes(y = 0), position = position_dodge(0.9)) +
  ylab("Loading coefficient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(scales = "free", space = "free_x")

ggsave(paste0(save_plots,"PLSDA_loadings.pdf"), plot=PLSDA_loadings)
saveRDS(PLSDA_loadings, paste0(save_plots,"PLSDA_loadings.rds"))


##################################################################
##                     Basic gPLS-DA model                      ##
##################################################################

#List of biomarkers in order of groups (based on Fran's grouping): 
#first 8 are liver, the next 10 are metabolic, next 2 immune, next 5 endocrine, final 3 kidney. 
#X_cuts defines these cuts. 
groups_fran = c('Alanine.aminotransferase','Alkaline.phosphatase','Aspartate.aminotransferase',
                'Direct.bilirubin','Gamma.glutamyltransferase','Total.bilirubin','Total.protein',
                'Albumin','Apolipoprotein.A','Apolipoprotein.B','Cholesterol','Glucose',
                'Glycated.haemoglobin.HbA1c.','HDL.cholesterol','LDL.direct','Lipoprotein.A',
                'Triglycerides','Urate','C.reactive.protein','IGF.1','Calcium','Phosphate',
                'SHBG','Testosterone','Vitamin.D','Creatinine','Cystatin.C','Urea')

X_fran = X[, groups_fran]
X_cuts_fran = c(8, 18, 20, 25)

#keepX is number of groups to keep  
gPLSDA <- gPLSda(X_fran, y, ncomp = 1, ind.block.x = X_cuts_fran, keepX = 1)
gPLSDA$loadings$X[gPLSDA$loadings$X != 0, ]



##################################################################
##             Fitting calibrated sPLS-DA model                 ##
##################################################################

#Create training and test sets 
smp_size <- floor(0.8*nrow(bio.cov))

#set the seed to same as in calibration to get the same train/test split
set.seed(123)
train_ind <- sample(seq_len(nrow(bio.cov)), size = smp_size)

train <- bio.cov[train_ind, ]
test <- bio.cov[-train_ind, ]

#Select all biomarkers from bio.cov for X
X_train = train[, 2:29]
X_test = test[, 2:29]
y_train = train$CVD_status
y_test = test$CVD_status


#Sparse plsda model; keepX is number of parameters to keep (9 had lowest misclassification rate in calibration). 
#The second line returns the variables selected.
sPLSDA <- splsda(X_train, y_train, ncomp=1, mode='regression', keepX=9)
sPLSDA$loadings$X[sPLSDA$loadings$X != 0, ]

#Predicting on the test set 
y_pred <- predict(sPLSDA, newdata = X_test)
fitted = y_pred$class$max.dist
table(fitted)
  

##################################################################
##                      Stability analyses                      ##
##################################################################

#Creating a heatmap of the selection of variables, over 100 iterations each selecting a different training/test set 
source("pls_functions.R")
set.seed(1)
Stability_results = StabilityPlot(X = X, Y = y, NIter = 100)

#pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE,
 #        display_numbers = TRUE, 
  #       filename = paste0(save_plots,"PLS_stability.pdf"),
   #      height = 5, width = 10)

#Calculating proportion of times each variable was selected
PropSelected = colSums(Stability_results)/28
PropSelected = as.data.frame(PropSelected)
library(dplyr)
PropSelected = tibble::rownames_to_column(PropSelected, "Biomarker")

stab_plot = PropSelected %>% 
  ggplot(aes(x = reorder(Biomarker, PropSelected))) +
  geom_linerange(aes(ymin = 0, ymax = PropSelected)) +
  coord_flip() + xlab("Biomarker") + ylab("Proportion selected")

ggsave(paste0(save_plots,"Stability_plot.pdf"), plot=stab_plot)
saveRDS(stab_plot, paste0(save_plots,"Stability_plot.rds"))



##################################################################
##             Fitting calibrated sgPLS-DA model                ##
##################################################################

#keepX is number of groups to keep; alpha is sparsity parameter (calibration not run yet)
sgPLSDA <- sgPLSda(X_fran, y, ncomp = 1, ind.block.x = X_cuts_fran, keepX = 3, alpha.x = 0.9)
sgPLSDA$loadings$X
sgPLSDA$loadings$X[sgPLSDA$loadings$X != 0, ]



#################################################################
##            Visualising the loadings coefficients            ##
#################################################################

#This plot visualises the loadings coefficients obtained from the fitted sPLSDA model

results = data.frame(cbind(Biomarker = colnames(X), Loadings = sPLSDA$loadings$X))

colnames(results)[2] = 'Loadings'
results$minLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) min(0, x)))
results$maxLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) max(0, x)))

sPLSDA_loadings = results %>% ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
                                        ymax = maxLoad))+
  geom_linerange(stat = "identity", position = position_dodge(0.9))+
  geom_point(aes(y = 0), position = position_dodge(0.9)) +
  ylab("Loading coefficient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(scales = "free", space = "free_x")

ggsave(paste0(save_plots,"sPLSDA_loadings.pdf"), plot=sPLSDA_loadings)
saveRDS(sPLSDA_loadings, paste0(save_plots,"sPLSDA_loadings.rds"))


#This plot visualises the loadings coefficients obtained from both sPLSDA and sgPLSDA models
#(will run this once I have calibrate the sgPLS)
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

sgPLSDA_loadings = results %>% ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
                                          ymax = maxLoad, color = Model)) +
  geom_linerange(stat = "identity", position = position_dodge(0.9)) +
  geom_point(aes(y = 0), position = position_dodge(0.9)) +
  ylab("Loading coefficients") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(cols = vars(belong_to), scales = "free", space = "free_x")

ggsave(paste0(save_plots,"sgPLSDA_loadings.pdf"), plot=sgPLSDA_loadings)
saveRDS(sgPLSDA_loadings, paste0(save_plots,"sgPLSDA_loadings.rds"))



#################################################################
##            Misclassification rates by subtype               ##
#################################################################

#Aim: copmute misclassification rates of PLS models by subtype of CVD.
#From full data, selected 5 most frequent CVD subtypes (based on ICD10 code for death):
#G459, I209, I219, I251, I639 (these all have over 600 cases). 

cvd_icd <- cov %>% select(ID, cvd_final_icd10)
bio.icd <- merge(bio, cvd_icd, by='ID')
y = bio.cov$CVD_status

#Create the stratified datasets 
for (subtype in c("G459","I209", "I219", "I249","I251", "I639")) {
  Xtemp = X[X$cvd_final_icd10 %in% subtype, ]
  Ytemp = y[X$cvd_final_icd10 %in% subtype]
  assign(paste0("X_", subtype), Xtemp)
  assign(paste0("Y_", subtype), Ytemp)
}

#Computing the misclassification rate for the calibrated sPLS-DA and sgPLS-DA models 
#And then creating a plot of these misclassification rates by CVD subtype
y_pred <- predict(sPLSDA, newdata = X)
fitted = y_pred$class$max.dist
table(fitted)

y_pred_g = predict(sgPLSDA, newdata = X)
fitted_g = y_pred_g$class$max.dist


#First, a plot of misclassification rates for sPLS-DA alone (as I have not calibrated the sgPLS-DA model)
mis_rate = data.frame(cbind(Prediction=fitted, subtype = bio.icd$cvd_final_icd10, 
                            truth=bio.cov$CVD_status))

levels(mis_rate$subtype) = c(levels(mis_rate$subtype),"Control")
mis_rate$subtype = replace(mis_rate$subtype, which(is.na(mis_rate$subtype)),
                           "Control")

mis_rate$IncorrectClass = !(mis_rate$comp1 == mis_rate$truth)

mis_rate_plot <- mis_rate %>% filter(subtype %in% c("Control","G459","I209", "I219", "I249",
                                                    "I251", "I639")) %>% 
  group_by(subtype) %>% 
  summarise(rate = mean(IncorrectClass)) %>% 
  arrange(subtype)

splsda_stratified <- mis_rate_plot %>% ggplot(aes(x = subtype, ymin = 0, ymax = rate)) + 
  geom_linerange(stat = "identity", position = position_dodge(0.9)) + 
  scale_color_brewer(palette = "Set1") + 
  ylab("Misclassification Rate")

ggsave(paste0(save_plots,"sPLSDA_stratified.pdf"), plot=splsda_stratified)
saveRDS(splsda_stratified, paste0(save_plots,"sPLSDA_stratified.rds"))


#Now for both sPLS-DA and sgPLS-DA 

mis_rate = data.frame(rbind(
  cbind(Prediction=fitted, model="sPLSDA", 
        subtype = bio.icd$cvd_final_icd10, truth=bio.cov$CVD_status),
  cbind(Prediction=fitted_g, model="sgPLSDA", 
        subtype = bio.icd$cvd_final_icd10, truth=bio.cov$CVD_status)))

levels(mis_rate$subtype) = c(levels(mis_rate$subtype),"Control")
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


