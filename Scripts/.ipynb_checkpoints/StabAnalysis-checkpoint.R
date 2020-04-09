
library(ggplot2)
library(dplyr)

platform = Sys.info()['sysname']
if (platform == 'Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

bio.imp = readRDS(paste0(data_folder,"bioImputedKNN.rds"))[,-1]

save_plots = paste0(save_plots,"PenalisedReg/ArrayJob/")

lasso.stab = readRDS(paste0(save_plots, "lassoStab1.rds"))
quietly(sapply(2:100, FUN =  function(x) {
  assign("lasso.stab", rbind(lasso.stab, readRDS(paste0(save_plots, "lassoStab",x,".rds"))),
         envir = .GlobalEnv)}))
lasso.stab = data.frame(lasso.stab)
colnames(lasso.stab) = colnames(bio.imp)

lasso.prop = apply(lasso.stab, 2, FUN = function(x) {
  sum(x != 0)/length(x)
}) 

lasso.mean = apply(lasso.stab, 2, mean)

lasso.sd = apply(lasso.stab, 2, sd)

resStabAnalysis = data.frame(
  Biomarker = colnames(bio.imp),
  PropSelected = lasso.prop,
    BetaMean = lasso.mean,
    BetaSD = lasso.sd
)

fig = resStabAnalysis %>% 
  ggplot(aes(x = reorder(Biomarker, PropSelected)))+
  geom_linerange(aes(ymin = 0, ymax = PropSelected))+
  coord_flip()+xlab("Biomarkers")+ylab("Proportion selected")
fig

ggsave(paste0(save_plots, "StabAnalysisLasso.pdf"))
saveRDS(fig, paste0(save_plots, "StabAnalysisLasso.pdf"))
