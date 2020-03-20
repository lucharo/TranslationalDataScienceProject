# LASSO, RIDGE AND E-NET ANALYSIS, BASED ON PRACTICAL 3

## cehck which variables were selected for lasso, enet


# NEED TO ADD CV/SUBSAMPLING KIND OF THING AND CHECK WHICH
# BIOS SELECTED ALL THE TIME

# LASSO AND PLS MAY BE MORE FOR EXPLORATORY ANALYSIS IN A WAY

library(glmnet)
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

ifelse(!dir.exists(file.path(save_plots, "PenalisedReg/")),
       dir.create(file.path(save_plots, "PenalisedReg/")), FALSE)
save_plots = paste0(save_plots,"PenalisedReg/")

bio = readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
bio.CVD = merge(bio,cov[,c("ID","CVD_status")],
                by="ID")
rownames(bio.CVD) = bio.CVD[,1]
bio.CVD = bio.CVD[,-1]
bio.CVD$CVD_status = as.factor(bio.CVD$CVD_status)

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

best.lam = "lambda.1se"

model_selector = function(alpha, plot = TRUE, save = T){
  set.seed(100) 
  # run cross validation
  model.cv <- cv.glmnet(X.train, y.train, alpha = alpha,
                     family = "binomial")
  
  betas = coef(model.cv, s = best.lam)[-1]
  ## Plots:
  #' 1. Calibration plot -> Deviance vs Lambda
  #' 2. Beta value for each covariate
  if (plot){
    if (save){
      # 1. Calibration plot
      pdf(file = paste0(save_plots,"CalibrationPlotAlpha",
                        as.character(alpha), ".pdf"), onefile = F)
      plot(model.cv)
      dev.off()
      # 2. Beta coefficients
      ggplot(data = data.frame(BetaValue = betas,
                               Biomarker = colnames(X)),
             aes(x = Biomarker, y = BetaValue))+
        geom_bar(stat = "identity")+coord_flip()+theme_minimal()
      ggsave(paste0(save_plots, "BetaValuesSingleRunAlpha",
                    as.character(alpha), ".pdf"))
    }
    plot(model.cv)
    
    print(ggplot(data = data.frame(BetaValue = betas,
                             Biomarker = colnames(X)),
           aes(x = Biomarker, y = BetaValue))+
      geom_bar(stat = "identity")+coord_flip()+theme_minimal())
  }

  
  bestlam = model.cv[[best.lam]]
  
  model_pred = predict(model.cv, s = bestlam,
                       newx = X.test )
  testMSE = mean((model_pred - y.test)^2)
  
  results = t(data.frame(c(round(alpha,2), round(model.cv$lambda.min,2),
                           round(model.cv$lambda.1se,2),
                           # get number of coefs with best
                           sum(coef(model.cv, s = bestlam)!=0)-1,
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
print(total_results)


