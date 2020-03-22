

library(ggplot2)
library(dplyr)
library(tidyr)

cluster = 0
platform = Sys.info()['sysname']
if (platform=='Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

save_plots = paste0(save_plots,"BHS/")

Mantej = readRDS(paste0(save_plots, "ScoresMantej.rds"))
Paper = readRDS(paste0(save_plots, "ScoresPaper.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
mini.cov = cov[,c("BS2_all", "ID", "CVD_status")]

myBHS = merge(Mantej, Paper, by = "ID")
colnames(myBHS)[2:3] = c("Mantej", "Paper")

allBHS = merge(myBHS, mini.cov, by = "ID")
# all(!duplicated(Mantej$ID))

allBHS %>% gather(key = "BHS", value = "value", -c(ID, CVD_status)) %>%
  ggplot(aes(x = as.factor(CVD_status), y = value)) +
  geom_boxplot()+facet_wrap(~BHS)
                  
                  
                  
                  
                  
                  
                  
