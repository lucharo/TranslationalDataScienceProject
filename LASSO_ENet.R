# LASSO, RIDGE AND E-NET ANALYSIS, BASED ON PRACTICAL 3

library(glmnet)

bio = readRDS("data/preprocessed/bioImputed.rds")
cov = readRDS("data/preprocessed/covProcessed.rds")
bio.CVD = merge(bio,cov[,"CVD_status"],
                by="row.names",all.x=TRUE)
rownames(bio.CVD) = bio.CVD$Row.names
bio.CVD = bio.CVD[,-1]

y = select(bio.CVD, CVD_status)
X = select(bio.CVD, -CVD_status)

set.seed(3141592) 
train.index = sample(1:nrow(X), 0.8*nrow(X))
X.train = as.matrix(X[train.index,])
y.train = as.numeric(y[train.index,1])
X.test = as.matrix(X[-train.index,])
y.test = as.numeric(y[-train.index,1])


###########################################################################
###########################################################################
###                                                                     ###
###                    RIDGE, LASSO AND ELASTIC NETS                    ###
###                                                                     ###
###########################################################################
###########################################################################

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
  model = cv.glmnet(x = X.train, y = y.train,
                    alpha = alpha, family = "binomial")
  with(model, cvm[which.min(lambda - lambda.1se)]) }

alpha.opt = optimise(cvm, c(0, 1))


## END
total_results = rbind(model_selector(0),
                      model_selector(1),
                      model_selector(alpha.opt$minimum))
total_results


