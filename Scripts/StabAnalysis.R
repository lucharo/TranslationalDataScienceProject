
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

bio.imp = readRDS(paste0(data_folder,"bioImputedKNN.rds"))[,-1]
save_plots = paste0(save_plots, "PenalisedReg/")
read_data = paste0(save_plots,"ArrayJob/")


lasso.stab = readRDS(paste0(read_data, "lassoStab1.rds"))
quietly(sapply(2:100, FUN =  function(x) {
  assign("lasso.stab", rbind(lasso.stab, readRDS(paste0(read_data, "lassoStab",x,".rds"))),
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
  ggplot(aes(x = reorder(Biomarker, PropSelected), color = as.factor(sign(BetaMean))))+
  geom_linerange(aes(ymin = 0, ymax = PropSelected))+
  # geom_col(aes(y = PropSelected))+
  coord_flip()+xlab("Biomarkers")+ylab("Proportion selected")
fig

ggsave(paste0(save_plots, "StabAnalysisLasso.pdf"))
saveRDS(fig, paste0(save_plots, "StabAnalysisLasso.rds"))

fig = resStabAnalysis %>% filter(PropSelected>0.5) %>% 
  ggplot(aes(x = reorder(Biomarker, BetaMean)))+
  geom_linerange(aes(ymin = 0, ymax = PropSelected,
                     color = as.factor(sign(BetaMean))), show.legend = F)+
  # geom_col(aes(y = PropSelected))+
  coord_flip()+xlab("Biomarkers")+ylab("Proportion selected")+
  theme_minimal()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+scale_color_brewer(palette = "Set1")
fig

fig2 = resStabAnalysis %>% filter(PropSelected>0.5)%>%
  ggplot(aes(x = reorder(Biomarker, BetaMean), exp(BetaMean)))+geom_point()+
  geom_linerange(aes(ymin = exp(BetaMean-BetaSD), ymax = exp(BetaMean+BetaSD)))+
  # geom_col(aes(y = PropSelected))+
  # coord_flip()+
  xlab("Biomarkers")+ylab("Odds Ratio")+coord_flip()+theme_minimal()
fig2

inOne = ggarrange(fig2, fig, widths = c(2,1), ncol = 2, common.legend = T)

inOne = annotate_figure(inOne,
                top = text_grob("Odd ratios of CVD incidence for biomarkers selected\n more than 50% of the time after 100 perturbation cycles"))

ggsave(paste0(save_plots, "ORStabAnalysis.pdf"))
saveRDS(fig, paste0(save_plots, "ORStabAnalysis.rds"))






