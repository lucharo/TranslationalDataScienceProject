# Stability analysis LASSO

library(glmnet)
library(ggplot2)
library(dplyr)
library(parallel)

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

############# Take args ################################
args = commandArgs(trailingOnly = TRUE)

seed = as.numeric(args[1])
print(seed)
#########################################################

print("Setting up data...")
t0 = Sys.time()
bio = readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
bio.CVD = merge(bio,cov[,c("ID","CVD_status")],
                by="ID")
rownames(bio.CVD) = bio.CVD[,1]
bio.CVD = bio.CVD[,-1]
bio.CVD$CVD_status = as.factor(bio.CVD$CVD_status)

y = dplyr::select(bio.CVD, CVD_status)
X = dplyr::select(bio.CVD, -CVD_status)

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
  ## add test score and return 
  return(coef.sub)          
}


################ 1 NODE SEVERAL CORES ###############
# niter = 100 
# print("Setting up cluster...")
# cl = makeCluster(detectCores()-1, type="FORK")
# print("Cluter set up")
# lasso.stab = parSapply(cl = cl, 1:niter, FUN = LassoSub, Xdata = as.matrix(X),
#                     Ydata = as.matrix(y))

# stopCluster(cl)
# print("Stability analysis finished.")

######## SEVERAL NODES FEW CORES ################

lasso.stab = LassoSub(k = seed, Xdata = as.matrix(X), Ydata = as.matrix(y))

ifelse(!dir.exists(file.path(save_plots, "ArrayJob/")),
       dir.create(file.path(save_plots, "ArrayJob/")), FALSE)
save_plots = paste0(save_plots,"ArrayJob/")

saveRDS(lasso.stab, paste0(save_plots, "lassoStab",as.character(seed),".rds"))







