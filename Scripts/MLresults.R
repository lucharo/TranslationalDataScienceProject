# Analyse ML results

library(ggplot2)

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

save_plots = paste0(save_plots, "MLAnalysis/")
read_data = paste0(save_plots,"ArrayJob/")

all.files = list.files(read_data)
results = readRDS(paste0(read_data, all.files[1]))

quietly(lapply(all.files[2:length(all.files)],
               FUN =  function(x) {
                 print(x)
                 assign("results",
                        rbind(results,
                              readRDS(paste0(read_data, x))),
                        envir = .GlobalEnv)}))

results %>% filter(data != "logBio") %>%
ggplot(aes(x = model, y = reorder(data, auc),
           fill = auc ,label = round(auc,4)))+
  geom_tile()+
  geom_text()+scale_fill_distiller(direction = 1)


