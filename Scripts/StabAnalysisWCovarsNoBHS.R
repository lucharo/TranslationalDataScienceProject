
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

save_plots = paste0(save_plots, "PenalisedRegWCov/")
read_data = paste0(save_plots,"ArrayJob/")


lasso.stab = readRDS(paste0(read_data, "lassoStabCovs1.rds"))
coefs = lasso.stab[[1]]
aucs = lasso.stab[[2]]
quietly(sapply(2:100, FUN =  function(x) {
  file = readRDS(paste0(read_data, "lassoStabCovs",x,".rds"))
  assign("coefs", rbind(coefs, file[[1]]),
         envir = .GlobalEnv)
  print(file[[2]])
  if (is.numeric(file[[2]])){
    assign("aucs", rbind(aucs, file[[2]]),
           envir = .GlobalEnv)
  } else { # sometimes results not stored as auc values but as ROC objects, weird
    assign("aucs", rbind(aucs, file[[2]]@y.values[[1]]),
           envir = .GlobalEnv)
  }

  }))

lasso.stab = data.frame(coefs)

lasso.prop = apply(lasso.stab, 2, FUN = function(x) {
  sum(x != 0)/length(x)
}) 

lasso.mean = apply(lasso.stab, 2, mean)

lasso.sd = apply(lasso.stab, 2, sd)

auc.mean = mean(aucs)
auc.sd = sd(aucs)
resStabAnalysis = data.frame(
  Variable = colnames(lasso.stab),
  PropSelected = lasso.prop,
    BetaMean = lasso.mean,
    BetaSD = lasso.sd
)

forTable = resStabAnalysis
forTable$Variable = str_replace_all(forTable$Variable, "\\.", " ")
forTable$Variable = str_replace_all(forTable$Variable, "_", " ")
forTable$Variable = str_replace_all(forTable$Variable, "Glycated haemoglobin HbA1c ", "HbA1c")
forTable$OR = round(exp(forTable$BetaMean),3)
forTable$LowerCI = round(exp(forTable$BetaMean-forTable$BetaSD),3)
forTable$UpperCI = round(exp(forTable$BetaMean+forTable$BetaSD),3)

forTable = forTable %>% select(-c(BetaMean, BetaSD))

forTable = forTable %>% arrange(-OR, -PropSelected)
forTable

resStabAnalysis$Variable = str_replace_all(resStabAnalysis$Variable, "\\.", " ")
resStabAnalysis$Variable = str_replace_all(resStabAnalysis$Variable, "Glycated haemoglobin HbA1c", "HbA1c")
resStabAnalysis$Variable = str_replace_all(resStabAnalysis$Variable, "Alanine aminotransferase", "Alanine\naminotransferase")
resStabAnalysis$Variable = str_replace_all(resStabAnalysis$Variable, "Alkaline phosphatase", "Alkaline\nphosphatase")
replace = c("# of medicines","no_medicines", "Age recruitment" , "age_recruitment 0 0",
            "BMI" , "BMI_5cl_2", "Male" , "gender Male", "Never smoker" , "smok_ever_2 no",
            "Education level" , "qual2", "Female", "gender Female")
for (index in seq(1,length(replace),2)){
  resStabAnalysis$Variable = str_replace(resStabAnalysis$Variable, replace[index+1], replace[index])
}


fig = resStabAnalysis %>% 
  ggplot(aes(x = reorder(Variable, PropSelected), color = as.factor(sign(BetaMean))))+
  geom_linerange(aes(ymin = 0, ymax = PropSelected))+
  # geom_col(aes(y = PropSelected))+
  coord_flip()+xlab("Variables")+ylab("Proportion selected")+
  theme(text = element_text(size = 15))
fig

ggsave(paste0(save_plots, "StabAnalysisLasso.pdf"))
saveRDS(fig, paste0(save_plots, "StabAnalysisLasso.rds"))

fig = resStabAnalysis %>% filter(PropSelected>0.5) %>% 
  ggplot(aes(x = reorder(Variable, BetaMean)))+
  geom_linerange(aes(ymin = 0, ymax = PropSelected,
                     color = as.factor(sign(BetaMean))), show.legend = F, size = 4)+
  # geom_col(aes(y = PropSelected))+
  coord_flip()+xlab("Variables")+ylab("Prop. selected")+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+scale_color_brewer(palette = "Set1")+
  theme(text = element_text(size = 20))+scale_y_continuous(labels = c("0","", "0.5", "", "1"))

fig2 = resStabAnalysis %>% filter(PropSelected>0.5)%>%
  ggplot(aes(x = reorder(Variable, BetaMean), exp(BetaMean)))+geom_point(size = 4)+
  geom_linerange(aes(ymin = exp(BetaMean-BetaSD), ymax = exp(BetaMean+BetaSD)), size = 2)+
  # geom_col(aes(y = PropSelected))+
  # coord_flip()+
  xlab("Variables")+ylab("Odds Ratio")+coord_flip()+theme_bw()+
  theme(text = element_text(size = 20))

inOne = ggarrange(fig2, fig, widths = c(3,1), ncol = 2, common.legend = T)
inOne
# inOne = annotate_figure(inOne,
#                 top = text_grob("Odd ratios of CVD incidence for Variables selected\n more than 50% of the time after 100 perturbation cycles"))

ggsave(paste0(save_plots, "ORStabAnalysis.pdf"))
saveRDS(fig, paste0(save_plots, "ORStabAnalysis.rds"))






