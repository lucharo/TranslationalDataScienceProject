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
suppressPackageStartupMessages(library(stringr))


cluster = 1

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = results_folder = "../FULLResults/"
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

bio.cov <- cov %>% 
  dplyr::select(ID, CVD_status) %>%
  merge(bio, by='ID')


##################################################################
##                     Basic PLS-DA model                       ##
##################################################################

#Select all biomarkers from bio.cov for X
X = bio.cov[, 3:30]
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
X_train = train[, 3:30]
X_test = test[, 3:30]
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

#No cases predicted
  

##################################################################
##                      Stability analyses                      ##
##################################################################

#Creating a heatmap of the selection of variables
#100 iterations each selecting a different training/test set 
source("pls_functions.R")
set.seed(1)
Stability_results = StabilityPlot(X = X, Y = y, NIter = 100)

saveRDS(Stability_results, paste0(save_plots,"Stability_results.rds"))

pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, 
         filename = paste0(save_plots,"PLS_stability.pdf"),
         height = 7, width = 10)


#Calculating proportion of times each variable was selected
PropSelected = colSums(Stability_results)/28
PropSelected = as.data.frame(PropSelected)
PropSelected = tibble::rownames_to_column(PropSelected, "Biomarker")

Prop50 = filter(PropSelected, PropSelected > 0.5)

stab_plot = Prop50 %>% 
  ggplot(aes(x = reorder(Biomarker, PropSelected))) +
  geom_linerange(aes(ymin = 0, ymax = PropSelected)) +
  #scale_x_continuous(limits = c(0.50, 1.00)) +
  coord_flip() + xlab("Biomarker") + ylab("Proportion selected")

ggsave(paste0(save_plots,"Stability_plot.pdf"), plot=stab_plot)
saveRDS(stab_plot, paste0(save_plots,"Stability_plot.rds"))



##################################################################
##             Fitting calibrated sgPLS-DA model                ##
##################################################################

#keepX is number of groups to keep; alpha is sparsity parameter 
#Calibration gave optimum number of groups = 2 and alpha = 0.5
#(but no cases were predicted)
sgPLSDA <- sgPLSda(X_fran, y, ncomp = 1, ind.block.x = X_cuts_fran, 
                   keepX = 2, alpha.x = 0.5)
sgPLSDA$loadings
#Since the Y2 is negative, the loadings should be reversed in the plot later
sgPLSDA$loadings$X[sgPLSDA$loadings$X != 0, ]



#################################################################
##            Visualising the loadings coefficients            ##
#################################################################

#calibrated sPLS on full data (not on training set as above) for plotting loadings 
sPLSDA <- splsda(X, y, ncomp=1, mode='regression', keepX=9)


#This plot visualises the loadings coefficients obtained from the fitted sPLSDA model

results = data.frame(cbind(Biomarker = colnames(X), Loadings = sPLSDA$loadings$X))

results = results %>%
  mutate(belong_to = ifelse(Biomarker %in% groups_fran[1:8], "Liver",
                            ifelse(Biomarker %in% groups_fran[9:18], "Metabolic",
                                   ifelse(Biomarker %in% groups_fran[19:20], "Immune",
                                          ifelse(Biomarker %in% groups_fran[21:25], "Endocrine",
                                                 "Kidney")))))

colnames(results)[2] = 'Loadings'
results$minLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) min(0, x)))
results$maxLoad = as.numeric(sapply(as.vector(results$Loadings), function(x) max(0, x)))

sPLSDA_loadings = results %>% ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
                                        ymax = maxLoad))+
  geom_linerange(stat = "identity", position = position_dodge(0.9)) +
  geom_point(aes(y = 0), position = position_dodge(0.9)) +
  ylab("Loading coefficient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(rows = vars(belong_to), scales = "free", space = "free_y") +
  theme(strip.text.y = element_text(angle = 0)) +
  coord_flip()

ggsave(paste0(save_plots,"sPLSDA_loadings.pdf"), plot=sPLSDA_loadings, height = 6)
saveRDS(sPLSDA_loadings, paste0(save_plots,"sPLSDA_loadings.rds"))


#This plot visualises the loadings coefficients obtained from both sPLSDA and sgPLSDA models
results_both = data.frame(rbind(
  cbind(Biomarker = colnames(X),
        Model = 'sPLSDA',
        Loadings = sPLSDA$loadings$X),
  cbind(Biomarker = colnames(X_fran),
        Model = 'sgPLSDA',
        Loadings = -(sgPLSDA$loadings$X))
))

results_both = results_both %>%
  mutate(belong_to = ifelse(Biomarker %in% groups_fran[1:8], "Liver",
                            ifelse(Biomarker %in% groups_fran[9:18], "Metabolic",
                                   ifelse(Biomarker %in% groups_fran[19:20], "Immune",
                                          ifelse(Biomarker %in% groups_fran[21:25], "Endocrine",
                                                 "Kidney")))))

colnames(results_both)[3] = 'Loadings'
results_both$minLoad = as.numeric(sapply(as.vector(results_both$Loadings), 
                                         function(x) min(0, x)))
results_both$maxLoad = as.numeric(sapply(as.vector(results_both$Loadings), 
                                         function(x) max(0, x)))

results_both$Biomarker = str_replace_all(results_both$Biomarker, "\\.", " ")
results_both$Biomarker = str_replace_all(results_both$Biomarker, 
                                         "Glycated haemoglobin HbA1c", "HbA1c")
results_both$Biomarker = str_replace_all(results_both$Biomarker, 
                                         "Gamma glutamyltransferase", "GGT")
results_both$Biomarker = str_replace_all(results_both$Biomarker, 
                                         "Alanine aminotransferase", "ALT")
results_both$Biomarker = str_replace_all(results_both$Biomarker, 
                                         "Aspartate aminotransferase", "AST")
results_both$Biomarker = str_replace_all(results_both$Biomarker, 
                                           "Alkaline phosphatase", "ALP")

sgPLSDA_loadings = results_both %>% ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
                                            ymax = maxLoad, color = Model)) +
    geom_linerange(stat = "identity", position = position_dodge(0.9)) +
    geom_point(aes(y = 0), position = position_dodge(0.9)) +
    ylab("Loading coefficients") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_brewer(palette = "Set1") +
    facet_grid(rows = vars(belong_to), scales = "free", space = "free_y") +
    theme(strip.text.y = element_text(angle = 0)) +
    coord_flip()

ggsave(paste0(save_plots,"sgPLSDA_loadings.pdf"), plot=sgPLSDA_loadings, height = 7.5)
saveRDS(sgPLSDA_loadings, paste0(save_plots,"sgPLSDA_loadings.rds"))



#################################################################
##            Misclassification rates by subtype               ##
#################################################################

#Aim: copmute misclassification rates of PLS models by subtype of CVD.
#From full data, selected 5 most frequent CVD subtypes (based on ICD10 code for death):
#G459, I209, I219, I251, I639 (these all have over 600 cases). 

bio.icd <- cov %>% 
  dplyr::select(ID, CVD_status, cvd_final_icd10) %>%
  merge(bio, by='ID')

#Computing the misclassification rate for the calibrated sPLS-DA and sgPLS-DA models 
#And then creating a plot of these misclassification rates by CVD subtype
y_pred <- predict(sPLSDA, newdata = X)
fitted = y_pred$class$max.dist
table(fitted)


#First, a plot of misclassification rates for sPLS-DA alone (as I have not calibrated the sgPLS-DA model)
mis_rate = data.frame(cbind(Prediction=fitted, subtype = bio.icd$cvd_final_icd10, 
                            truth=bio.cov$CVD_status))

levels(mis_rate$subtype) = c(levels(mis_rate$subtype),"Control")
mis_rate$subtype = replace(mis_rate$subtype, which(is.na(mis_rate$subtype)),
                           "Control")

mis_rate$IncorrectClass = !(mis_rate$comp1 == mis_rate$truth)

mis_rate_plot <- mis_rate %>% filter(subtype %in% c("Control","G459","I209","I219",
                                                    "I251","I639")) %>% 
  group_by(subtype) %>% 
  summarise(rate = mean(IncorrectClass)) %>% 
  arrange(subtype)

splsda_stratified <- mis_rate_plot %>% ggplot(aes(x = subtype, ymin = 0, ymax = rate)) + 
  geom_linerange(stat = "identity", position = position_dodge(0.9)) + 
  scale_color_brewer(palette = "Set1") + 
  ylab("Misclassification Rate")

ggsave(paste0(save_plots,"sPLS_mislcassification.pdf"), plot=splsda_stratified)
saveRDS(splsda_stratified, paste0(save_plots,"sPLS_mislcassification.rds"))




#################################################################
##                     Stratified analyses                     ##
#################################################################

X_strat = bio.icd[, 4:31]
y_strat = bio.icd$CVD_status

#Create the stratified datasets 
for (subtype in c("G459","I209","I219","I251","I639")) {
  Xtemp = X_strat[bio.icd$cvd_final_icd10 %in% subtype | is.na(bio.icd$cvd_final_icd10), ]
  Ytemp = y_strat[bio.icd$cvd_final_icd10 %in% subtype | is.na(bio.icd$cvd_final_icd10)]
  assign(paste0("X_", subtype), Xtemp)
  assign(paste0("Y_", subtype), Ytemp)
}

#Make sPLS models for each subtype
for (subtype in c("G459","I209", "I219","I251","I639")) {
  X = eval(parse(text = paste0("X_", subtype)))
  Y = eval(parse(text = paste0("Y_", subtype)))
  sPLSDA_strat <- splsda(X, Y, keepX = 9, ncomp = 1, mode = "regression")
  assign(paste0("sPLSDA_", subtype), sPLSDA_strat)
}


#Plot loadings 
results = data.frame(rbind(
  cbind(Biomarker = colnames(X),
        Subtype = 'G459',
        Loadings = sPLSDA_G459$loadings$X),
  cbind(Biomarker = colnames(X),
        Subtype = 'I209',
        Loadings = sPLSDA_I209$loadings$X),
  cbind(Biomarker = colnames(X),
        Subtype = 'I219',
        Loadings = sPLSDA_I219$loadings$X),
  cbind(Biomarker = colnames(X),
        Subtype = 'I251',
        Loadings = sPLSDA_I251$loadings$X),
  cbind(Biomarker = colnames(X),
        Subtype = 'I639',
        Loadings = sPLSDA_I639$loadings$X)
))


results_strat = results %>%
  mutate(belong_to = ifelse(Biomarker %in% groups_fran[1:8], "Liver",
                            ifelse(Biomarker %in% groups_fran[9:18], "Metabolic",
                                   ifelse(Biomarker %in% groups_fran[19:20], "Immune",
                                          ifelse(Biomarker %in% groups_fran[21:25], "Endocrine",
                                                 "Kidney")))))

colnames(results_strat)[3] = 'Loadings'
results_strat$minLoad = as.numeric(sapply(as.vector(results_strat$Loadings), 
                                   function(x) min(0, x)))
results_strat$maxLoad = as.numeric(sapply(as.vector(results_strat$Loadings), 
                                          function(x) max(0, x)))

results_strat$Biomarker = str_replace_all(results_strat$Biomarker, "\\.", " ")
results_strat$Biomarker = str_replace_all(results_strat$Biomarker, 
                                         "Glycated haemoglobin HbA1c", "HbA1c")
results_strat2$Biomarker = str_replace_all(results_strat2$Biomarker, 
                                           "Alkaline phosphatase", "ALP")

strat_loadings = results_strat %>% ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
                                          ymax = maxLoad, color = Subtype)) +
  geom_linerange(stat = "identity", position = position_dodge(0.9)) +
  geom_point(aes(y = 0), position = position_dodge(0.9)) +
  ylab("Loading coefficients") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(cols = vars(belong_to), scales = "free", space = "free_x")

ggsave(paste0(save_plots,"sPLSDA_stratified.pdf"), plot=strat_loadings)
saveRDS(strat_loadings, paste0(save_plots,"sPLSDA_stratified.rds"))



##Redo this plot excluding 0 coefficients for simplicity 

results_strat2 = results %>%
  mutate(belong_to = ifelse(Biomarker %in% groups_fran[1:8], "Liver",
                            ifelse(Biomarker %in% groups_fran[9:18], "Metabolic",
                                   ifelse(Biomarker %in% groups_fran[19:20], "Immune",
                                          ifelse(Biomarker %in% groups_fran[21:25], "Endocrine",
                                                 "Kidney")))))

colnames(results_strat2)[3] = 'Loadings'
results_strat2$Loadings = as.character(results_strat2$Loadings)
results_strat2$Loadings = as.numeric(results_strat2$Loadings)

biomarkers_non0 = results_strat2 %>% 
  group_by(Biomarker) %>% 
  summarise(SumCoefs = sum(abs(Loadings))) %>%
  filter(SumCoefs != 0)
  
idx_keep = results_strat2$Biomarker %in% biomarkers_non0$Biomarker
results_strat2 = results_strat2[idx_keep, ]

results_strat2$minLoad = as.numeric(sapply(as.vector(results_strat2$Loadings), 
                                          function(x) min(0, x)))
results_strat2$maxLoad = as.numeric(sapply(as.vector(results_strat2$Loadings), 
                                          function(x) max(0, x)))

results_strat2$Biomarker = str_replace_all(results_strat2$Biomarker, "\\.", " ")
results_strat2$Biomarker = str_replace_all(results_strat2$Biomarker, 
                                          "Glycated haemoglobin HbA1c", "HbA1c")
results_strat2$Biomarker = str_replace_all(results_strat2$Biomarker, 
                                           "Gamma glutamyltransferase", "GGT")
results_strat2$Biomarker = str_replace_all(results_strat2$Biomarker, 
                                           "Alkaline phosphatase", "ALP")

strat_loadings2 = results_strat2 %>% 
  ggplot(aes(x = Biomarker, y = 0, ymin = minLoad,
             ymax = maxLoad, color = Subtype)) +
  geom_linerange(stat = "identity", position = position_dodge(0.9)) +
  geom_point(aes(y = 0), position = position_dodge(0.9)) +
  ylab("Loading coefficients") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(rows = vars(belong_to), scales = "free", space = "free_y") +
  theme(strip.text.y = element_text(angle = 0)) +
  coord_flip()

ggsave(paste0(save_plots,"sPLSDA_strat_non0.pdf"), plot=strat_loadings2, height=6.5)
saveRDS(strat_loadings2, paste0(save_plots,"sPLSDA_strat_non0.rds"))
