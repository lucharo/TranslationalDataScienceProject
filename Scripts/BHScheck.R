

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

library(ROCR)


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
colnames(mini.cov)[colnames(mini.cov) == "BS2_all"]="BS2"

myBHS = merge(Mantej, Paper, by = "ID")
myBHS = merge(myBHS, Barbara, by = "ID")
colnames(myBHS)[2:4] = c("B", "A", "C")

allBHS = merge(myBHS, mini.cov, by = "ID")
# all(!duplicated(Mantej$ID))

allBHS %>% gather(key = "BHS", value = "value", -c(ID, CVD_status)) %>%
  ggplot(aes(x = as.factor(CVD_status), y = value)) +
  geom_boxplot(aes(fill = BHS), show.legend = F)+
  facet_grid(~BHS)+
  stat_compare_means(method = "t.test", paired = F, label = "p.signif",
                     label.x = 1.5, label.y = 1.1)+
  stat_summary(geom = "point", shape = 23, fun.data = 'mean_se')+
  xlab("CVD status")+ylab("BHS score")+
  # ggtitle("BHS score by CVD status")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  theme(text = element_text(size = 24))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))

################# BHS evaluation #################################

# NEED TO consider sampling methods like under/oversampling and potentially bootstrapping
# to always take into account the cases and account for the big imbalance

y = as.numeric(allBHS$CVD_status)
X = allBHS[,-c(1, ncol(allBHS))]

k.folds = 5
## set up for CV
folds <- cut(seq(1,nrow(X)),breaks=k.folds,labels=FALSE)

best_auc = 0
best_BHS = 0

######### CV
BHSes = colnames(X)

for (BHS in 1:4){
  auc.list = c()
  for (i in 1:k.folds){
    valIndexes <- which(folds==i)
    X.val <- X[valIndexes, BHS]
    X.train <- X[-valIndexes, BHS]
    
    y.val = y[valIndexes]
    y.train = y[-valIndexes]
    
    train = data.frame(y.train = y.train, X.train = X.train)
    val = data.frame(X.val = X.val)
    colnames(val) = "X.train"
    
    mod = glm(y.train~X.train, data = train, family = "binomial")
    prdct = round(predict.glm(mod, newdata =val, type = "response"))
    pred.objct = prediction(prdct, y.val)
    auc = performance(pred.objct, "auc")
    auc.list = append(auc.list, auc@y.values[[1]])
  }
  mean.auc = mean(auc.list)
  if (mean.auc>best_auc){
    best_auc = mean.auc
    best_BHS= BHSes[BHS]
  }
  print(mean.auc)
}

print(best_BHS)

                  
                  
                  
                  
                  
