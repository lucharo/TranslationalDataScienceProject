# EDA COV

library(tidyverse)
library(naniar)

cluster = 0
platform = Sys.info()['sysname']
if (platform == 'Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

ifelse(!dir.exists(file.path(save_plots, "EDA/")), dir.create(file.path(save_plots, "EDA/")), FALSE)
save_plots = paste0(save_plots,"EDA/")

save.results = function(figure, name, ggsv = T){
  if (ggsv){
    ggsave(paste0(save_plots,name,".pdf"))
  }
  saveRDS(figure, paste0(save_plots,name,".rds"))
}

cov = readRDS(paste0(data_folder,"covProcessed.rds"))
cov.original = readRDS(paste0("../FULLDATA/Covariates_full.rds"))

cov <- cov %>% 
  select(-c("vit_status","dc_cvd_st","age_cl", "stop","stop_cvd",
            "age_CVD", "cvd_final_icd10", "primary_cause_death_ICD10",
            "cvd_death", "cvd_incident", "cvd_prevalent"))
cov.original = cov.original %>% select(
  colnames(cov.original)[colnames(cov.original) %in% colnames(cov)])

fig = cov.original %>%
  dplyr::group_by(CVD_status) %>%
  miss_var_summary() %>%
  ggplot(aes(CVD_status,
             variable,
             fill = pct_miss)) +
  geom_tile() +
  viridis::scale_fill_viridis(name = "% Miss") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle("Missing data patterns by CVD status for all biomarkers")
fig



cov %>% ggplot(aes(x = CVD_status, y = age_recruitment.0.0, fill = gender))+geom_boxplot()+
  stat_compare_means(aes(group = gender), label = "p.format", method = "t.test")+theme_minimal()+
  scale_fill_brewer(palette = "Set1")+xlab("CVD status")+ylab("Age at recruitment")









