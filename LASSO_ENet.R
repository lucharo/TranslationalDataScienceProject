# LASSO, RIDGE AND E-NET ANALYSIS, BASED ON PRACTICAL 3

## cehck which variables were selected for lasso, enet


library(glmnet)

bio = readRDS("data/preprocessed/bioImputed.rds")
cov = readRDS("data/preprocessed/covProcessed.rds")
bio.CVD = merge(bio,cov[,"CVD_status"],
                by="row.names",all.x=TRUE)
rownames(bio.CVD) = bio.CVD$Row.names
bio.CVD = bio.CVD[,-1]

y = dplyr::select(bio.CVD, CVD_status)
X = dplyr::select(bio.CVD, -CVD_status)

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

best.lam = "lambda.min"

model_selector = function(alpha){
  set.seed(100) 
  # run cross validation
  model <- cv.glmnet(X.train, y.train, alpha = alpha,
                     family = "binomial")
  
  plot(model)
  
  bestlam = model[[best.lam]]
  
  model_pred = predict(model, s = bestlam,
                       newx = X.test )
  testMSE = mean((model_pred - y.test)^2)
  
  results = t(data.frame(c(round(alpha,2), round(model$lambda.min,2),
                           round(model$lambda.1se,2),
                           # get number of coefs with best
                           sum(coef(model, s = bestlam)!=0)-1,
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
  # pick the error (given by cvm) with the best lambda (lamda.1se or lambda.min)
  # in theory we'd want lambda == lamba.min or lambda.1se but for some reason 
  # lambda does not containg any value equal to lambda.min or lambda.1se, 
  # therefore we do the maths:
  #     lambda == lambda.min
  # ==> lambda - lambda.min == 0 
  # this difference could be negative but we're interested in the lambda 
  # closer to lambda min (smaller absolute difference)
  # ==> lambda for which abs(lambda-lambda.min) is minimum (closer to 0)
  with(model, cvm[which.min(abs(lambda - lambda.min))]) }

alpha.opt = optimise(cvm, c(0, 1))


## END
total_results = rbind(model_selector(0),
                      model_selector(1),
                      model_selector(alpha.opt$minimum))
total_results


