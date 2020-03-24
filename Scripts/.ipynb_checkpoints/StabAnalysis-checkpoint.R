
library(ggplot2)
library(dplyr)

lasso.prop = apply(lasso.stab, 1, FUN = function(x) {
  sum(x != 0)/length(x)
}) 
names(lasso.prop) = colnames(X)

resStabAnalysis = data.frame(
  Biomarker = colnames(X),
  PropSelected = lasso.prop
)

fig = resStabAnalysis %>% 
  ggplot(aes(x = reorder(Biomarker, PropSelected)))+
  geom_linerange(aes(ymin = 0, ymax = PropSelected))+
  coord_flip()+xlab("Biomarkers")+ylab("Proportion selected")
fig

ggsave(paste0(save_plots, "StabAnalysisLasso.pdf"))
saveRDS(fig, paste0(save_plots, "StabAnalysisLasso.pdf"))
