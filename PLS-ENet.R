## PLS/Elastic Net attempt on TDS data


library(glmnet)
library(lme4)
# Aim of this script is to replicate the work from practical 3 and 4 on the
# TDS dataset

# Maybe also do univariate analysis using mixed models effect (maybe for
# duplicate measurements)

###########################################################################
###########################################################################
###                                                                     ###
###                    RIDGE, LASSO AND ELASTIC NETS                    ###
###                                                                     ###
###########################################################################
###########################################################################

# load biomarkers dataset, ideally one of the already preprocessed ones
# so you dont have to run preprocessing again
biomarkers = readRDS("data/MARbiomarkers.rds")
y = select(biomarkers, CVD_status)
X = select(biomarkers, -CVD_status)

set.seed(3141592) 
train.index = sample(1:nrow(X), 0.8*nrow(X))
X.train = as.matrix(X[train.index,])
y.train = as.numeric(y[train.index,1])
X.test = as.matrix(X[-train.index,])
y.test = as.numeric(y[-train.index,1])

model_selector = function(alpha){
  set.seed(100) 
  # run cross validation
  model <- cv.glmnet(X.train, y.train, alpha = alpha,
                     family = "binomial")
  
  plot(model)
  
  bestlam = model$lambda.1se
  
  model_pred = predict(model, s = bestlam,
                       newx = X.test )
  testMSE = mean((model_pred - y.test)^2)
  
  results = t(data.frame(c(round(alpha,2), round(model$lambda.min,2),
                           round(model$lambda.1se,2),
                           length(which(coef(model)!=0))-1,
                           testMSE)))
  colnames(results) = c("Alpha","lambda.min","lambda.1se",
                        "# of predictors", "Test MSE")
  rownames(results) = ifelse(alpha == 0, "Ridge",
                             ifelse(alpha == 1, "Lasso", "ElasticNet"))
  
  return(results)
}


# alpha optimisation for Elastic Net
cvm = function(alpha) {
  model = cv.glmnet(x = X.train, y = y.train, alpha = alpha)
  with(model, cvm[which.min(lambda - lambda.1se)]) }
alpha.opt = optimise(cvm, c(0, 1))


## END
total_results = rbind(model_selector(0),
                      model_selector(1),
                      model_selector(alpha.opt$minimum))
total_results



############################################################################
############################################################################
###                                                                      ###
###                     PARTIAL LEAST SQUARES MODELS                     ###
###                                                                      ###
############################################################################
############################################################################

