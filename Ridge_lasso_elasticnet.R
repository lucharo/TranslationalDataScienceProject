
# Load packages -----------------------------------------------------------

library(glmnet)
library(tidyverse)
library(caret)


# Training and Testing Set ------------------------------------------------

bio_imp = readRDS("data/preprocessed/bioImputed.rds")
cov_pro = readRDS("data/preprocessed/covProcessed.rds")

bio.CVD = merge(bio_imp,cov_pro[,"CVD_status"],
                by="row.names",all.x=TRUE)
rownames(bio.CVD) = bio.CVD$Row.names
bio.CVD = bio.CVD[,-1]

covariates_processed$cvd_final_icd10 <- NULL
covariates_processed$primary_cause_death_ICD10 <- NULL

covariates_processed$CVD_status <- ifelse(covariates_processed$CVD_status == "1", 1, 0)
Y <- as.numeric(unlist(covariates_processed[, "CVD_status"]))
X <- model.matrix(CVD_status~., covariates_processed)[,-1]

set.seed(100)
train <- sample(1:nrow(covariates_processed), 0.7 * nrow(covariates_processed))
test <- seq(1, nrow(covariates_processed))[-train]


# Ridge Regression --------------------------------------------------------

set.seed(100)
model_ridge <- cv.glmnet(X[train, ], Y[train], alpha =0, family  = "binomial")
plot(model_ridge)

model_ridge$lambda.min
model_ridge$lambda.1se
min(model_ridge$cvm)

id_min <- which(model_ridge$lambda == model_ridge$lambda.min)
model_ridge$cvm[id_min]

id_1se_ridge <- which(model_ridge$lambda == model_ridge$lambda.1se)
model_ridge$cvm[id_1se_ridge]

best1am_ridge <- model_ridge$lambda.1se

model_ridge_pred = predict(model_ridge,  s = best1am_ridge, newx = X[test,])
mean((model_ridge_pred - Y[test])^2)


# Lasso -------------------------------------------------------------------

set.seed((100))
model_lasso <- cv.glmnet(X[train,], Y[train], alpha = 1, family = "binomial")
plot(model_lasso)

model_lasso$lambda.min
model_lasso$lambda.1se
min(model_lasso$cvm)

id_1se_lasso <- which(model_lasso$lambda == model_lasso$lambda.1se)
model_lasso$cvm[id_1se_lasso]

bestlam_lasso <- model_lasso$lambda.1se
length(which(coef(model_lasso, s = bestlam_lasso) !=0)) -1

model_lasso_pred <- predict(model_lasso, s = bestlam_lasso, newx = X[test, ])
mean((model_lasso_pred - Y[test])^2)


# Elastic Net -------------------------------------------------------------

cvm = function(alpha) {
  model = cv.glmnet(x = X[train,], y = Y[train], alpha = alpha)
  with(model, cvm[which.min(lambda - lambda.1se)])
}
set.seed(100)
alpha.opt = optimize(cvm, c(0,1))
alpha.opt

model_enet = cv.glmnet(x = X, y = Y, alpha = alpha.opt$minimum, family = "binomial")  
plot(model_enet)

model_enet$lambda.min
model_enet$lambda.1se
min(model_enet$cvm)

id_1se_enet <- which(model_enet$lambda == model_enet$lambda.1se)
model_enet$cvm[id_1se_enet]

bestlam_enet <- model_enet$lambda.1se
length(which(coef(model_enet, s = bestlam_enet) !=0)) -1

model_enet_pred <- predict(model_enet, s = bestlam_enet, newx = X[test, ])
mean((model_enet_pred - Y[test])^2)
