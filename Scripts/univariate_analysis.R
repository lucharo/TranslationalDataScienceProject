library(glmnet)
library(lme4)
library(dplyr)
library(ggplot2)
# Aim of this script is to replicate the work from practical 3 and 4 on the
# TDS dataset

cluster = 1
platform = Sys.info()['sysname']
if (platform == 'Linux'){
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

bio.imp_cov = merge(bio.imp,cov,by="ID")
bio.imp_cov = bio.imp_cov[,-1]

bio_cov = merge(bio,cov,by="ID")
bio_cov = bio_cov[,-1]

confounders = c("age_CVD",
                # "BS2_all", excludng BS2 as it depends on the biomarker,
                # not meet requirements to be a valid confounder
                "qual2","smok_ever_2",
                "physical_activity_2", "alcohol_2", "BMI_5cl_2",
                "no_cmrbt_cl2", "no_medicines", "gender")


### ordered factors
tidyFactors = function(df){
  df$qual2 = factor(df$qual2, levels = c("low", "intermediate", "high"), ordered = T)
  df$alcohol_2 = factor(df$alcohol_2, levels = c("Non-drinker", "Social drinker",
                                                           "Moderate drinker", "Daily drinker"), ordered = T)
  df$BMI_5cl_2 = factor(df$BMI_5cl_2, levels = c("[18.5,25[", "[25,30[",
                                                           "[30,40[", ">=40"), ordered = T)
  df$no_cmrbt_cl2 = as.numeric(df$no_cmrbt_cl2)
  df$no_medicines = as.numeric(df$no_medicines)
  df$smok_ever_2 = factor(df$smok_ever_2, levels=c("no","yes"))
  df$physical_activity_2 = factor(df$physical_activity_2, levels=c("no","yes"))
  df$CVD_status = factor(df$CVD_status, levels = c("0","1"))
  df$gender = factor(df$gender, levels=c("Female","Male"))
  
  return(df)
}

bio_cov = tidyFactors(bio_cov)
bio.imp_cov = tidyFactors(bio.imp_cov)


bio_cov = bio_cov[, c("CVD_status", confounders, colnames(bio)[-1])]
bio.imp_cov = bio.imp_cov[, c("CVD_status", confounders, colnames(bio.imp)[-1])]


DoYouMatter = function(biomarker, data = bio.imp_cov){
  # I guess I would ideally not have to precise a data argument and allow
  # the data to be the global proteins.covars object, though lmer throws an
  # error
  
  
  formula = as.formula(paste("CVD_status",
                             paste(c(confounders, biomarker),
                                   collapse = "+"), sep = "~"))
  model = glm(formula, data = data,
              na.action = na.omit, family = "binomial")
  
  # first element is beta coefficient (log OR), second is p-value
  list(summary(model)$coefficients[biomarker, 1],
       summary(model)$coefficients[biomarker, 4],
        confint(model)[biomarker, 1:2])
}

Univariate.analysis = function(merged.dataset){
  
  
  bio.names = colnames(bio.imp)[-1]
  
  pvals = data.frame(
    matrix(unlist(lapply(bio.names,
                         function(x) 
                           DoYouMatter(x,data = merged.dataset))),
           ncol = 4, byrow = T),
    stringsAsFactors = F)
  
  Univ_biomarkers = data.frame("Biomarkers" = bio.names, pvals)
  
  colnames(Univ_biomarkers)[2:5] = c("OR", "p.value", "lower CI", "upper CI")
  
  Univ_biomarkers
}

univ_nonimputed = Univariate.analysis(bio_cov)
univ_imputed = Univariate.analysis(bio.imp_cov)

##################################################################
##                           Plotting                           ##
##################################################################

fortable = rbind(cbind(univ_nonimputed, Data = "Not-imputed"),
                 cbind(univ_imputed, Data = "KNN")) %>% arrange(p.value)
fortable$Biomarkers = str_replace_all(fortable$Biomarkers, "\\.", " ")
fortable$Biomarkers = str_replace_all(fortable$Biomarkers, "Glycated haemoglobin HbA1c", "HbA1c")
fortable$OR = exp(fortable$OR)
fortable$`lower CI` = exp(fortable$`lower CI`)
fortable$`upper CI` = exp(fortable$`upper CI`)

readr::write_csv(fortable, paste0(save_plots,"UnivariateAnalysis.csv"))

data = rbind(cbind(univ_nonimputed, Data = "Not-imputed"),
             cbind(univ_imputed, Data = "KNN")) %>% arrange(p.value) %>% head(30)

data$Biomarkers = str_replace_all(data$Biomarkers, "\\.", " ")
data$Biomarkers = str_replace_all(data$Biomarkers, "Glycated haemoglobin HbA1c", "HbA1c")
data$Biomarkers = str_replace_all(data$Biomarkers, "Alanine aminotransferase", "Alanine\naminotransferase")
data$Biomarkers = str_replace_all(data$Biomarkers, "Alkaline phosphatase", "Alkaline\nphosphatase")


figure = data %>%
  ggplot(aes(x = reorder(Biomarkers, -log10(p.value)),
             y = -log10(p.value),
             shape = as.factor(
               ifelse(OR>0,"Positive",
                      "Negative")),
             #size = as.numeric(abs(OR)),
             color = Data,
             size = exp(OR)))+
  geom_point(position = position_dodge(0.9))+
  geom_hline(yintercept = -log10(0.05/length(colnames(bio))),
             linetype = "dashed")+
  # geom_text(aes(x = Biomarkers, y = -log10(p.value)+0.4,
  #               label = ifelse(Data == "MICE" &
  #                                -log10(p.value)> -log10(0.05/length(colnames(bio))),
  #                              colnames(bio)[Biomarkers], "")),
  #           angle = 90, show.legend = F, color = "black", size = 3)+
  scale_shape_manual(values = c(6,2))+
  scale_color_brewer(palette = "Set1")+
  labs(shape = "Sign of the\nassociation")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
    # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank()
    legend.position = c(0.8, 0.4),
    legend.background = element_rect(fill="white", 
                                     size=0.2, linetype="solid"),
    text = element_text(size = 13)
    )+
  ylim(0, NA)+
  ylab("-log10 of p-values")+xlab("Biomarkers")+
  labs(size = "Odds Ratio")+
  # theme(legend.position = "top")+
  ggtitle("Univariate analysis results for biomarkers")+
  coord_flip()
figure
##################################################################
##                         Saving Plots                         ##
##################################################################

ggsave(paste0(save_plots,"UnivariateAnalysis.pdf")) # as pdf
saveRDS(figure, paste0(save_plots,"UnivariateAnalysis.rds")) # as object



## sorting data by Odd ratio insteadof pvalue

# 
# figure2 = rbind(cbind(univ_nonimputed, Data = "Not-imputed"),
#                cbind(univ_imputed, Data = "KNN")) %>% arrange(desc(OR)) %>% head(20) %>%
#   ggplot(aes(x = reorder(Biomarkers, -exp(OR)),
#              y = -log10(p.value),
#              shape = as.factor(
#                ifelse(OR>0,"Positive",
#                       "Negative")),
#              #size = as.numeric(abs(OR)),
#              color = Data,
#              size = exp(OR)))+
#   geom_point(position = position_dodge(0.9))+
#   geom_hline(yintercept = -log10(0.05/length(colnames(bio))),
#              linetype = "dashed")+
#   # geom_text(aes(x = Biomarkers, y = -log10(p.value)+0.4,
#   #               label = ifelse(Data == "MICE" &
#   #                                -log10(p.value)> -log10(0.05/length(colnames(bio))),
#   #                              colnames(bio)[Biomarkers], "")),
#   #           angle = 90, show.legend = F, color = "black", size = 3)+
#   scale_shape_manual(values = c(6,2))+
#   scale_color_brewer(palette = "Set1")+
#   labs(shape = "Sign of the\nassociation")+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1),
#         # axis.text.x = element_blank(),
#         # axis.ticks.x = element_blank()
#   )+
#   ylim(0, NA)+
#   ylab("-log10 of p-values")+xlab("Biomarkers")+
#   labs(size = "Odds Ratio")+theme(legend.position = "top")+
#   ggtitle("Univariate analysis results for imputed and non imputed data")
# figure2
