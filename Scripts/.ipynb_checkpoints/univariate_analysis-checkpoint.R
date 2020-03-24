library(glmnet)
library(lme4)
library(dplyr)
library(ggplot2)
# Aim of this script is to replicate the work from practical 3 and 4 on the
# TDS dataset

cluster = 1

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

###########################################################################
###########################################################################
###                                                                     ###
###                UNIVARIATE ANALYSIS OF THE BIOMARKERS                ###
###                                                                     ###
###########################################################################
###########################################################################

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))
bio.imp = readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))

bio.imp_cov = merge(bio.imp,cov,by="row.names",all.x=TRUE)
rownames(bio.imp_cov) = bio.imp_cov$Row.names
bio.imp_cov = bio.imp_cov[,-1]

bio_cov = merge(bio,cov,by="row.names",all.x=TRUE)
rownames(bio_cov) = bio_cov$Row.names
bio_cov = bio_cov[,-1]

confounders = c("age_CVD","BS2_all", "qual2","smok_ever_2",
                "physical_activity_2", "alcohol_2", "BMI_5cl_2",
                "no_cmrbt_cl2", "no_medicines")

DoYouMatter = function(biomarker, data = bio.imp_cov){
  # I guess I would ideally not have to precise a data argument and allow
  # the data to be the global proteins.covars object, though lmer throws an
  # error
  
  formula = as.formula(paste("CVD_status",
                             paste(c(confounders, biomarker),
                                   collapse = "+"), sep = "~"))
  model = glm(formula, data = data,
              na.action = na.omit, family = "binomial")
  
  # first element is beta coefficient (OR), second is p-value
  list(summary(model)$coefficients[biomarker, 1],
       summary(model)$coefficients[biomarker, 4])
}

Univariate.analysis = function(merged.dataset){
  
  
  bio.names = colnames(bio.imp)
  
  pvals = data.frame(
    matrix(unlist(lapply(bio.names,
                         function(x) 
                           DoYouMatter(x,data = merged.dataset))),
           ncol = 2, byrow = T),
    stringsAsFactors = F)
  
  Univ_biomarkers = data.frame("Biomarkers" = bio.names, pvals)
  
  colnames(Univ_biomarkers)[2:3] = c("OR", "p.value")
  
  Univ_biomarkers
}

univ_nonimputed = Univariate.analysis(bio_cov)
univ_imputed = Univariate.analysis(bio.imp_cov)

##################################################################
##                           Plotting                           ##
##################################################################

figure = rbind(cbind(univ_nonimputed, Data = "Not-imputed"),
      cbind(univ_imputed, Data = "MICE")) %>% 
  ggplot(aes(x = Biomarkers,
             y = -log10(p.value),
             shape = as.factor(
               ifelse(OR>0,"Positive",
                      "Negative")),
             #size = as.numeric(abs(OR)),
             color = Data))+
  geom_point(position = position_dodge(0.9))+
  geom_hline(yintercept = -log10(0.05/length(colnames(bio))),
             linetype = "dashed")+
  scale_shape_manual(values = c(6,2))+
  scale_color_brewer(palette = "Set1")+
  labs(shape = "Sign of the\nassociation")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##################################################################
##                         Saving Plots                         ##
##################################################################

ggsave(paste0(save_plots,"UnivariateAnalysis.pdf")) # as pdf
saveRDS(figure, paste0(save_plots,"UnivariateAnalysis.rds")) # as object

