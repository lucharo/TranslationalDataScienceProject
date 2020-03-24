

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)


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
Barbara = readRDS(paste0(save_plots, "ScoresBarbara.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
mini.cov = cov[,c("BS2_all", "ID", "CVD_status")]

myBHS = merge(Mantej, Paper, by = "ID")
myBHS = merge(myBHS, Barbara, by = "ID")
colnames(myBHS)[2:4] = c("Mantej", "Paper", "Barbara")

allBHS = merge(myBHS, mini.cov, by = "ID")
# all(!duplicated(Mantej$ID))

allBHS %>% gather(key = "BHS", value = "value", -c(ID, CVD_status)) %>%
  ggplot(aes(x = as.factor(CVD_status), y = value)) +
  geom_boxplot()+
  facet_wrap(~BHS)+
  stat_compare_means(method = "t.test", paired = F,
                     label.x = 1.5, label.y = 0.9,
                     comparisons = list(c(0, 1)))+
  stat_summary(geom = "point", shape = 23, fun.data = 'mean_se')
                  
                  
                  
                  
                  
                  
                  
