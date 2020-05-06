# Aim of this script is to calculate a polygenic risk score for each person

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(jtools)
library(scales)
library(tidyr)
library(ggsignif)
library(ROCR)


cluster = 1

if (cluster == 1){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = save_data = "../FULLResults/PRS/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

snp = readRDS(paste0(data_folder,"snpProcessed.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
snp.info = readRDS(paste0(data_folder,"snpInfo.rds"))


##################################################################
##                        Computing PRS                         ##
##################################################################

#Added this for cluster as snp file is smaller in cluster than in toy? 
#162 snps vs 177 - since those with more than 90% similarity were removed?
snps = colnames(snp)
snp.info = snp.info[snp.info$markername %in% snps, ]

#Check that the snps are in the same order in snp and snp.info ([-1] to remove the ID column)
stopifnot(all(colnames(snp)[-1] == snp.info$markername))

#Extracting the beta values (which are in the correct order based on above check)
betas <- snp.info$beta

#Element-wise multiplication of each row of snp by the betas
#(margin=2 specifies it is rows of the matrix; each row is a person)
#This multiples each no. of snp copies by the beta coefficient for that snp 
PRS_matrix = sweep(snp[,-1], MARGIN=2, betas, `*`)


#Calculate the weighted sum for each person (na.rm=TRUE to ignore missing values)
PRS_sums = rowSums(PRS_matrix, na.rm=TRUE)

#Calculate the average (divide by the number of SNPs with data available)
PRS = PRS_sums/rowSums(!is.na(snp))

PRS = cbind(snp$ID, PRS)
PRS = as.data.frame(PRS)
colnames(PRS) = c("ID", "PRS")
saveRDS(PRS, paste0(save_data,"PolygenicRiskScore.rds"))


#################################################################
##                Comparison between CVD status                ##
#################################################################

## Merge cov and PRS
cov.prs <- merge(cov, PRS, by="ID")

#Normalise the PRS
cov.prs$PRS = rescale(cov.prs$PRS, to = c(-1, 1)) 

#Density plot
den_plot <- ggplot(cov.prs) + 
  geom_density(aes(x=PRS, color=CVD_status, fill=CVD_status), alpha=0.2, size=0.25) +
  labs(x = 'Polygenic Risk Score', y = 'Density') + 
  scale_fill_discrete(name = "CVD status", labels = c("Controls", "Cases")) +
  guides(color = FALSE) + 
  theme_minimal()

ggsave(paste0(save_plots,"PRS_density.pdf"), den_plot)

#Boxplot of PRS by CVD status
prs_boxplot <- ggplot(cov.prs, aes(x=CVD_status, y=PRS)) +
  geom_boxplot() +
  theme_minimal() + 
  labs(x = 'CVD Status', y = 'Polygenic Risk Score')
ggsave(paste0(save_plots,"PRS_boxplot.pdf"), prs_boxplot)

#t-test - sig difference in mean PRS between groups
t_test = t.test(PRS ~ CVD_status, data=cov.prs)
saveRDS(t_test, paste0(save_plots,"PRS_ttest.rds"))



###Association between PRS and CVD - simple logistic regression model 

cov.prs$CVD_status <- as.numeric(cov.prs$CVD_status)
glm <- glm(CVD_status ~ as.numeric(PRS), data=cov.prs, family=binomial)
summary(glm)
odds <- round(cbind("odds" = exp(coef(glm)), exp(confint(glm))), 3)
saveRDS(odds, paste0(save_plots,"PRS_Odds.rds"))

glm_adj <- glm(CVD_status ~ PRS + gender + age_recruitment.0.0, data=cov.prs, family=binomial)
summary(glm_adj)
round(cbind("odds" = exp(coef(glm_adj)), exp(confint(glm_adj))), 3)


#Visualising the relationship between PRS and risk of CVD
pdf(paste0(save_plots,"PRS_EffectPlot.pdf"))
glm_plot <- effect_plot(glm, pred = PRS, data=cov.prs, interval = TRUE, int.width = 0.95) + 
  labs(x = 'Polygenic Risk Score', y = 'Probability of having CVD')
glm_plot
dev.off()
saveRDS(glm_plot, paste0(save_plots,"PRS_EffectPlot.rds"))


###Predicting CVD by PRS - 5-fold cross-validation logistic models -> AUC

X = cov.prs$PRS
y = cov.prs$CVD_status

k.folds = 5
folds <- cut(seq(1,nrow(X)), breaks=k.folds, labels=FALSE)

best_auc = 0
best_PRS = -1

######### CV
auc.list = c()

for (i in 1:k.folds){
  valIndexes <- which(folds == i)
  X.val <- X[valIndexes]
  X.train <- X[-valIndexes]
    
  y.val = y[valIndexes]
  y.train = y[-valIndexes]
    
  train = data.frame(y.train = y.train, X.train = X.train)
  val = data.frame(X.val = X.val)
  colnames(val) = "X.train"
    
  mod = glm(y.train ~ X.train, data = train, family = "binomial")
  pred = round(predict.glm(mod, newdata = val, type = "response"))
  pred.object = prediction(pred, y.val)
  auc = performance(pred.object, "auc")
  auc.list = append(auc.list, auc@y.values[[1]])
}

mean.auc = mean(auc.list)


