# opti analysis

library(ggplot2)
library(dplyr)
library(ggpubr)

platform = Sys.info()['sysname']
if (platform == 'Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

save_plots = paste0(save_plots, "knnOpti/")
read_data = paste0(save_plots,"ArrayJob/")


# load data
list.files(read_data)
opti.results = readRDS(paste0(read_data, "knnOpti1.rds"))
absolute = opti.results[[1]]
perParam = opti.results[[2]]
perParam$k = seq(1,20,2)
perParam = perParam %>% pivot_longer(-k)

quietly(sapply(2:50, FUN =  function(x) {
  file =  readRDS(paste0(read_data, "knnOpti",x,".rds"))
  temp.perParam = file[[2]]
  temp.perParam$k = seq(1,20,2)
  temp.perParam = temp.perParam %>% pivot_longer(-k)
  assign("absolute", rbind(absolute, file[[1]]),
         envir = .GlobalEnv)
  assign("perParam", rbind(perParam, temp.perParam),
         envir = .GlobalEnv)
  
}))


# Absolute errors

absolute.mean = colMeans(absolute)
absolute.sd = apply(absolute, 2, sd)

res.absol = data.frame(mean = absolute.mean, sd = absolute.sd,
                       k = seq(1,20,2))

res.absol %>% ggplot(aes(x = k, y = mean, ymax = mean+sd, ymin = mean-sd))+geom_pointrange()

# Per Param error
perParam$k = as.factor(perParam$k)
sumup = perParam %>% group_by(name, k) %>% dplyr::summarise(value = mean(value))
sumup %>%
  ggplot(aes(x = k, y = reorder(name, value), fill = value, label = round(value,2)))+geom_tile(color = "black")+
  geom_text()+scale_fill_distiller(direction = 1, palette = "Reds")+ylab("Biomarkers")+
  xlab("Amount of nearest neighbours considered(k)")+
  labs(fill = "Mean estimated\nRMSE")+
  theme(text = element_text(size = 12), axis.text.y = element_text(angle = 45))
