# Stability analysis LASSO

library(glmnet)
library(ggplot2)
library(caret)
library(dplyr)
library(parallel)
library(ROCR)

platform = Sys.info()['sysname']
if (platform == 'Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

PRSdf = readRDS(paste0(save_plots, "PRS/PolygenicRiskScore.rds"))

bestBHS = readRDS(paste0(save_plots, "BHS/ScoresPaper.rds"))

ifelse(!dir.exists(file.path(save_plots, "NoBioPenalisedRegWCov/")),
       dir.create(file.path(save_plots, "NoBioPenalisedRegWCov/")), FALSE)
save_plots = paste0(save_plots,"NoBioPenalisedRegWCov/")

############# Take args ################################
args = commandArgs(trailingOnly = TRUE)

seed = as.numeric(args[1])
if (length(args)==0) seed = 1
print(seed)
#########################################################

print("Setting up data...")
t0 = Sys.time()
cov = readRDS(paste0(data_folder,"covProcessed.rds"))

cov <- cov %>% # we will leave CVD_stattus in cov and let it be split in the
  # analysis function so that the sampling can be done and order is kept
  select(-c("vit_status","dc_cvd_st","age_cl", "stop","stop_cvd",
            "age_CVD", "cvd_final_icd10", "primary_cause_death_ICD10",
            "cvd_death", "cvd_incident", "cvd_prevalent", "BS2_all"))

##################################################################
##                       ONE HOT ENCODING                       ##
##################################################################

cov$smok_ever_2 = as.factor(cov$smok_ever_2)
cov$physical_activity_2 = as.factor(cov$physical_activity_2)
cov$gender = as.factor(cov$gender)

onehotvars = dummyVars("~.", data = cov[,c("smok_ever_2","physical_activity_2", "gender")])
onehotvars = data.frame(predict(onehotvars, newdata = cov[,c("smok_ever_2",
                                                             "physical_activity_2",
                                                             "gender")]))
# Delete one hot encoded variables and add the encoded versions
cov = cov %>% select(-c(smok_ever_2,physical_activity_2, gender)) %>% cbind(onehotvars)
# remove onehot vars for memory management
remove(onehotvars)

##################################################################
##                 Factors to ordinal variables                 ##
##################################################################

# factors to ordinal variables
cov$qual2 = factor(cov$qual2, levels = c("low", "intermediate", "high"), ordered = T)
cov$alcohol_2 = factor(cov$alcohol_2, levels = c("Non-drinker", "Social drinker",
                                                 "Moderate drinker", "Daily drinker"), ordered = T)
cov$BMI_5cl_2 = factor(cov$BMI_5cl_2, levels = c("[18.5,25[", "[25,30[",
                                                 "[30,40[", ">=40"), ordered = T)
cov$no_cmrbt_cl2 = as.numeric(cov$no_cmrbt_cl2)
cov$no_medicines = as.numeric(cov$no_medicines)



##################################################################
##                   Merging all data sources                   ##
##################################################################

data = merge(cov,bestBHS, by = "ID")

data$CVD_status = as.factor(data$CVD_status)

data = merge(data, PRSdf, by = "ID")

data$PRS = scales::rescale(data$PRS, to = c(-1,1))

rownames(data) = data$ID
data = data %>% select(-ID)

#################################################################
##                          X-y split                          ##
#################################################################

y = dplyr::select(data, CVD_status)
X = dplyr::select(data, -CVD_status)

print("Data setup")
print(Sys.time()-t0)

best.lam = "lambda.1se"

## Sensitivity analysis - for lasso
LassoSub = function(k = 1, Xdata, Ydata, family = "binomial",
                    penalty.factor = NULL){
  if (is.null(penalty.factor)){
    penalty.factor = rep(1, ncol(Xdata)) 
  }
  set.seed(k)
  s = sample(nrow(Xdata), size = 0.8 * nrow(Xdata))
  Xsub = Xdata[s, ]
  Ysub = Ydata[s, 1] 
  
  model.sub = cv.glmnet(x = Xsub, y = Ysub, alpha = 1,
                        family = family,
                        penalty.factor = penalty.factor)
  coef.sub = coef(model.sub, s = best.lam)[-1]
  
  ### Test
  Xtest = Xdata[-s,]
  ytest = Ydata[-s,1]
  best.prdct = predict(model.sub, s = "lambda.1se", newx = Xtest)
  pred.objct = ROCR::prediction(best.prdct, ytest)
  auc = ROCR::performance(pred.objct, "auc")
  names(coef.sub) = colnames(X)
  ## add test score and return 
  return(list(coef.sub,auc@y.values[[1]]))          
}


######## SEVERAL NODES FEW CORES ################

lasso.stab = LassoSub(k = seed, Xdata = data.matrix(X), Ydata = as.matrix(y))

ifelse(!dir.exists(file.path(save_plots, "ArrayJob/")),
       dir.create(file.path(save_plots, "ArrayJob/")), FALSE)
save_plots = paste0(save_plots,"ArrayJob/")

saveRDS(lasso.stab, paste0(save_plots, "NoBiolassoStabCovs",as.character(seed),".rds"))







