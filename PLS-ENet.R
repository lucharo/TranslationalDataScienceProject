## PLS/Elastic Net attempt on TDS data


library(glmnet)
library(lme4)
# Aim of this script is to replicate the work from practical 3 and 4 on the
# TDS dataset

# Maybe also do univariate analysis using mixed models effect (maybe for
# duplicate measurements)

# load biomarkers dataset, ideally one of the already preprocessed ones
# so you dont have to run preprocessing again
biomarkersCVD = readRDS("data/MARbiomarkers.rds")
y = select(biomarkersCVD, CVD_status)
X = select(biomarkersCVD, -CVD_status)

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

###########################################################################
###########################################################################
###                                                                     ###
###                UNIVARIATE ANALYSIS OF THE BIOMARKERS                ###
###                                                                     ###
###########################################################################
###########################################################################

bio.covs = readRDS("data/bio_and_covs.rds")

confounders = c("age_CVD","BS2_all", "qual2","smok_ever_2",
                "physical_activity_2", "alcohol_2", "BMI_5cl_2",
                "no_cmrbt_cl2", "no_medicines")

bio.names = colnames(select(biomarkersCVD, -CVD_status))

DoYouMatter = function(molecule, data = bio.covs){
  # I guess I would ideally not have to precise a data argument and allow
  # the data to be the global proteins.covars object, though lmer throws an
  # error
  
  formula0 = as.formula(paste(molecule,
                              paste(confounders, collapse = "+"), sep = "~"))
  m0 = lm(formula0, data = data, na.action = na.omit)
  
  formula1 = as.formula(paste(molecule,
                              paste(c(confounders, "CVD_status"),
                                    collapse = "+"),
                              sep = "~"))
  m1 = lm(formula1, data = data, na.action = na.omit)
  
  anova(m0,m1)$`Pr(>Chisq)`[2]
}

pvals = unlist(lapply(bio.names, DoYouMatter))
FatANOVA = data.frame("Biomarkers" = bio.names, "p-value" = pvals)

