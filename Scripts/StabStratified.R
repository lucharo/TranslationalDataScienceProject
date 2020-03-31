#Prepare libraries and data 

suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(sgPLS))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))


cluster = 1

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_results = results_folder = "../FULLResults/StabStratified/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_results = "../results/StabStratified/"
}

bio <- readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov <- readRDS(paste0(data_folder,"covProcessed.rds"))

bio.icd <- cov %>% dplyr::select(ID, CVD_status, cvd_final_icd10) %>%
  inner_join(bio, by='ID')

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


#Stability analyses for each model 
source("pls_functions.R")
set.seed(1)
for (subtype in c("G459","I209","I219","I251","I639")) {
  print(subtype)
  MyStab = NULL
  X = eval(parse(text = paste0("X_", subtype)))
  Y = eval(parse(text = paste0("Y_", subtype)))
  Stability_results = StabilityPlot(X = X, Y = Y, NIter = 100)
  saveRDS(Stability_results, paste0(save_results,"Stability_",subtype,".rds"))
}



#################################################################
##         Stability analyses of of stratified models          ##
#################################################################

#Read in selection heatmaps 
G459 <- readRDS(paste0(results_folder,"Stability_G459.rds"))
I209 <- readRDS(paste0(results_folder,"Stability_I209.rds"))
I219 <- readRDS(paste0(results_folder,"Stability_I219.rds"))
I251 <- readRDS(paste0(results_folder,"Stability_I251.rds"))
I639 <- readRDS(paste0(results_folder,"Stability_I639.rds"))


#Calculating proportion of times each variable was selected in each of the models
for (subtype in c("G459","I209","I219","I251","I639")) {
  data = eval(parse(text = subtype))
  proportions = colSums(data)/28
  proportions = as.data.frame(proportions)
  proportions = tibble::rownames_to_column(proportions, "Biomarker")
  assign(paste0("PropSelected_", subtype), proportions)
}


#Make into a dataframe for plotting
stab_strat = as.data.frame(rbind(
  cbind(Biomarker = PropSelected_G459$Biomarker,
        Subtype = 'G459',
        PropSelected = PropSelected_G459$proportions),
  cbind(Biomarker = PropSelected_I209$Biomarker,
        Subtype = 'I209',
        PropSelected = PropSelected_I209$proportions),
  cbind(Biomarker = PropSelected_I219$Biomarker,
        Subtype = 'I219',
        PropSelected = PropSelected_I219$proportions),
  cbind(Biomarker = PropSelected_I251$Biomarker,
        Subtype = 'I251',
        PropSelected = PropSelected_I251$proportions),
  cbind(Biomarker = PropSelected_I639$Biomarker,
        Subtype = 'I639',
        PropSelected = PropSelected_I639$proportions)
), stringsAsFactors = FALSE)

stab_strat$PropSelected = as.numeric(stab_strat$PropSelected)


#Concatenate these into one dataframe 
#First, rename columns so they are identifiable after joining
#colnames(PropSelected_I459)[2] = "PropSelected_I459"
#colnames(PropSelected_I209)[2] = "PropSelected_I209"
#colnames(PropSelected_I219)[2] = "PropSelected_I219"
#colnames(PropSelected_I251)[2] = "PropSelected_I251"
#colnames(PropSelected_I639)[2] = "PropSelected_I639"

#Now join them
#PropSelected = cbind(PropSelected_I209, PropSelected_I219, PropSelected_I251,
#                     PropSelected_I639)
#Drop extra biomarker columns
#PropSelected = PropSelected[-c(3,5,7)]


#Plot those selected above 50% of the time 
stab_strat50 = filter(stab_strat, PropSelected > 0.5)

#Or biomarkers with average selection above 50%?
#select = stab_strat %>% group_by(Biomarker) %>% 
 # summarize(mean = mean(PropSelected)) %>%
  #filter(mean > 0.5)

stab_plot = stab_strat %>%
  ggplot(aes(x = reorder(Biomarker, PropSelected),
             ymin = 0, ymax = PropSelected, color=Subtype)) +
  geom_linerange(stat = "identity", position = position_dodge(0.9)) +
  coord_flip() + xlab("Biomarker") + ylab("Proportion selected") +
  scale_color_brewer(palette = "Set1")

ggsave(paste0(save_results,"StabStratified_plot.pdf"), plot=stab_plot)
#saveRDS(stab_plot, paste0(save_results,"StabStratified_plot.rds"))
