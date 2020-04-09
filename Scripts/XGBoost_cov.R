
# Installing and Loading Packages -----------------------------------------
if (!require(xgboost)) install.packages("xgboost")
library(xgboost)
if (!require(readr)) install.packages("readr")
library(readr)
if (!require(stringr)) install.packagaes("stringr")
library(stringr)
if (!require(caret)) install.packages("caret")
library(caret)
if (!require(car)) install.packages("car")
library(car)

library(tidyverse)

# for auc calculation
library(ModelMetrics)

# for svm
library(e1071)

cluster = 1
platform = Sys.info()['sysname']
if (platform == 'Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

#Load Data
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
bio = readRDS(paste0(data_folder,"bioImputed.rds"))
snp = readRDS(paste0(data_folder,"snpImputed.rds"))

# Preparing our Data and selecting features -------------------------------

# One stumbling block when getting started with the xgboost package in R is
# that you can't just pass it a dataframe. The core xgboost function requires data to be a matrix.

#To prepare our data, we have a number of steps we need to complete:
#remove information about the target variable from the training data
#reduce the amount of redundant information
#convert categorical information to a numeric format
#Split dataset into training and testing subsets
#Convert the cleaned dataframe to a Dmatrix


# Remove information about the target variable from the training d --------

#First let's remove the columns that have information on our target variable

#Let's create a new vector with labels - convert the CVD_status factor to an
# integer class starting at 0, as the first class should be 0; picky requirement
cov$CVD_status = as.integer(cov$CVD_status)

cov <- cov %>% # we will leave CVD_stattus in cov and let it be split in the
  # analysis function so that the sampling can be done and order is kept
  select(-c("vit_status","dc_cvd_st","age_cl", "stop","stop_cvd",
            "age_CVD", "cvd_final_icd10", "primary_cause_death_ICD10",
            "cvd_death", "cvd_incident", "cvd_prevalent"))

# Reduce the amount of redundant information ------------------------------

#Finally, I want to remove all the non-numeric variables,
# since a matrix can only hold numeric variables - if you try to create a
# matrix from a dataframe with non-numeric variables in it, it will
# conver them all into NA's and give you some warning messages.

# Luis: we want to keep ordered categorical variables as integers but OneHotEncoding
# on non-ordered categorical variables (e.g. smok_ever, physical activity) only those
# 2 surprisingly

cov$smok_ever_2 = as.factor(cov$smok_ever_2)
cov$physical_activity_2 = as.factor(cov$physical_activity_2)

onehotvars = dummyVars("~.", data = cov[,c("smok_ever_2","physical_activity_2")])
onehotvars = data.frame(predict(onehotvars, newdata = cov[,c("smok_ever_2",
                                                             "physical_activity_2")]))
# Delete one hot encoded variables and add the encoded versions
cov = cov %>% select(-c(smok_ever_2,physical_activity_2)) %>% cbind(onehotvars)
# remove onehot vars for memory management
remove(onehotvars)

# Convert dataframes into matrix -------------------------------------------

# cov <- as.matrix(cov)
# bio = as.matrix(as.numeric(bio))

# preferred xgboost data type = xgb.Dmatrix

#################################################################
##                  Function would start here                  ##
#################################################################

# Analysis = function(model, data, outcome = 'CVD_status'){
#   #reproducibility across datasets
#
#   set.seed(1247)
#   tr.index = sample(1:nrow(data), size = 0.7*nrow(data))
#   X = dplyr::select(data, -outcome)
#   y = dplyr::select(data, outcome)
#
#   y.train = y[tr.index,1]
#   X.train = X[tr.index,]
#
#   y.test = y[-tr.index,1]
#   X.test = X[-tr.index,]
#
#   if (model == "xgboost"){
#     X.train = as.matrix(X.train)
#     X.test = as.matrix(X.test)
#     train = xgb.DMatrix(data = X.train, label= y.train)
#     test <- xgb.DMatrix(data = X.test, label = y.test)
#
#     # nrounds is basically number of trees in forest
#     # this is basically a grid search
#     best_auc = 0
#     best_eta = 0
#     best_depth = 0
#     best_nround = 0
#     for (nround in c(5,10,20)){
#       #Apparently one should only tweak nrounds realistically or maybe yes
#       # https://www.kaggle.com/c/otto-group-product-classification-challenge/discussion/13053
#       # https://stackoverflow.com/questions/35050846/xgboost-in-r-how-does-xgb-cv-pass-the-optimal-parameters-into-xgb-train
#       # many other things to optimise: https://rdrr.io/cran/xgboost/man/xgb.train.html
#       for (max_depth in c(3,5,10,15)){
#         for (eta in c(0.1, 0.5, 1)){
#
#           # if the nfold is too tiny this function may give an error because the dataset
#           # fed into the model only contains negative samples aka (controls)
#
#           model.cv = xgb.cv(data = train, nfold = 5, nrounds = nround, objective = "binary:logistic",
#                       metrics = list("auc"), eta = eta, max_depth = max_depth, verbose = F)
#           if (max(model.cv[["evaluation_log"]][["test_auc_mean"]]) > best_auc){
#             best_eta = eta
#             best_depth = max_depth
#             best_nround = nround
#             best_auc = max(model.cv[["evaluation_log"]][["test_auc_mean"]])
#           }
#         }
#       }
#     }
#
#     best.model = xgb.train(data = train, nrounds = best_nround, objective = "binary:logistic",
#                            eta = best_eta, max_depth = best_depth, eval_metric = "auc")
#     prdct = predict(best.model, newdata = test)
#     pred <- prediction(prdct, getinfo(test,'label'))
#     perf <- performance(pred,"tpr","fpr")
#     auc = performance(pred, "auc")
#     plot(perf, coloured = T)
#     return(auc@y.values[[1]])
#
#   }
#
#   else if (model == "svm"){
#
#     ## set up for CV
#     folds <- cut(seq(1,nrow(X.train)),breaks=5,labels=FALSE)
#
#     best_auc = 0
#     best_kernel = 0
#     best_cost = 0
#     for (kernel in c("linear","radial")){
#       for(cost in c(1,5,10,15)){
#         ######### CV
#         auc.list = c()
#         for (i in 1:5){
#           valIndexes <- which(folds==i)
#           ValData <- X.train[valIndexes, ]
#           nonValData <- X.train[-valIndexes, ]
#
#           ValOutcome = y.train[valIndexes]
#           nonValOutcome = y.train[-valIndexes]
#
#
#           # details for this syntax :
#           #  https://stackoverflow.com/questions/9028662/predict-maybe-im-not-understanding-it
#           md.svm = svm(nonValOutcome ~ . , data = data.frame(nonValOutcome, nonValData),
#                        kernel = kernel, cost = cost)
#           prdct = predict(md.svm, newdata =  ValData)
#           pred.objct = prediction(prdct, ValOutcome)
#           auc = performance(pred.objct, "auc")
#           auc.list = append(auc.list, auc@y.values[[1]])
#         }
#         mean.auc = mean(auc.list)
#         if (mean.auc>best_auc){
#           best_auc = mean.auc
#           best_kernel = kernel
#           best_cost = cost
#         }
#       }
#     }
#
#     bst.svm = svm(y.train~., data = data.frame(y.train, X.train),
#                  kernel = best_kernel, cost = best_cost)
#     best.prdct = predict(bst.svm, X.test)
#     pred.objct = prediction(best.prdct, y.test)
#     auc = performance(pred.objct, "auc")
#     plot(performance(pred.objct, "tpr", "fpr"))
#     return(auc@y.values[[1]])
#   }
#
#   else if (model == "logit"){
#
#   }
# }
#

Analysis("xgboost", cov)

Analysis("svm", cov)

cov.bio = merge(cov,bio, by = "row.names")
rownames(cov.bio) = cov.bio$Row.names
cov.bio = cov.bio[,-1]

Analysis("svm", cov.bio)

Analysis("xgboost", cov.bio)

cov.bio.snp = merge(cov.bio, snp, by = "row.names")
rownames(cov.bio.snp) = cov.bio.snp$Row.names
cov.bio.snp = cov.bio.snp[,-1]

Analysis("svm", cov.bio.snp)

Analysis("xgboost", cov.bio.snp)

# more to come: +BHS_score, +PRS, +both ...






# Split dataset into Training and Testing Sets ----------------------------

#Get the numb 70/30 training test split
tr.index.cov <- sample(1:nrow(cov), size = 0.7*nrow(cov))
tr.index.bio = sample(1:nrow(cov), size = 0.7*nrow(cov))
#Training Data
cov.train <- cov[tr.index, ]
y.train <- Outcome[tr.index]

#Testing Data
cov.test <- cov[-tr.index, ]
test_labels <- Outcome[-tr.index]


# Convert the cleaned dataframe to a DMATRIX ------------------------------

# put our testing & training data into two seperates Dmatrixs objects
dtrain <- xgb.DMatrix(data = train_data, label= train_labels)
dtest <- xgb.DMatrix(data = test_data, label= test_labels)


# Training our model ------------------------------------------------------

#In order to train our model, we need to give it some information to start with. It needs to know:

#What training data to use. In this case, we've already put our data in a dmatrix and can just pass it that.
#The number of training rounds. This just means the number of times we're going to improve our naive model by adding additional models. In the figure above going around the circle once = one boosting round.
#What the objective function is. The objective function you pick will depend on the task you have. Because we're going to try and predict something that's binary (either humans get sick or they don't), we're going to use “binary:logistic”, which is logistic regression for binary (two-class) classification. By default, xgboost will do linear regression.

# train a model using our training data
model <- xgboost(data = dtrain, # the data
                 nround = 2, # max number of boosting iterations
                 objective = "binary:logistic")

# generate predictions for our held-out testing data
pred <- predict(model, dtest)

# get & print the classification error
err <- mean(as.numeric(pred > 0.5) != test_labels)
print(paste("test-error=", err))


# Tuning our model --------------------------------------------------------

#By default, the max depth of trees in xgboost is 6 - let's set it to 3 instead

# train an xgboost model
model_tuned <- xgboost(data = dtrain, # the data
                       max.depth = 3, # the maximum depth of each decision tree
                       nround = 2, # max number of boosting iterations
                       objective = "binary:logistic")

# generate predictions for our held-out testing data
pred <- predict(model_tuned, dtest)

# get & print the classification error
err <- mean(as.numeric(pred > 0.5) != test_labels)
print(paste("test-error=", err))

#Hence we have not over-fitted


# Examining our model -----------------------------------------------------

# plot them features! what's contributing most to our model?
xgb.plot.multi.trees(feature_names = names(covaraites_matrix),
                     model = model, fill = TRUE)
#Doesn't seem to work for some reason???

model_tree <- xgb.dump(model_tuned, with_stats = T)
model_tree[1:10]

#Get information about how important each feature is
names <- dimnames(dtrain)[[2]]
importance_matrix <- xgb.importance(names, model = model_tuned)
xgb.plot.importance(importance_matrix[1:10,])


# Alternative Method ------------------------------------------------------

xgb_params <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = nc)

watchlist <- list(train = dtrain, test = dtest)

#Extreme Gradeint Boosting Model
bst_model <- xgb.train(params = xgb_params,
                       data = dtrain,
                       nrounds = 100, #Maximum number of iterations
                       watchlist = watchlist )
#train-mlogloss = This prints out what the error was on the training data
#test-mlogloss = This prints out what the error was on the testing data
#In the first iteration there is a higher error in the testing data than the training data

#Training and test error plot
bst_model
#Makes use of the evaluation_log results to make it into a plot
e <- data.frame(bst_model$evaluation_log)
plot(e$iter, e$train_mlogloss, col = "blue")
#Fairly high error rate at the beginning but as the iterations increase the error rate significantly does too
#Now let's add the test error to this plot
lines(e$iter, e$test_mlogloss, col = "red")

#Possible over-fitting - need to adjust the model perhaps?

#To see the minimum value of the test-error
min(e$test_mlogloss)
#To know which iteration this minimum value belonged to
e[e$test_mlogloss == 0.179392, ] #Hence the 9th iteration was the lowest

#Now to optimise our model:
bst_model_1 <- xgb.train(params = xgb_params,
                       data = dtrain,
                       nrounds = 100, #Maximum number of iterations
                       watchlist = watchlist,
                       eta = 0.2) #Default value for eta is 0.3 and can range from 0 to 1 - this is a learning rate; so if this value is low it means it is more robust to overfitting - hence can avoid overfitting by keeping this value low)
bst_model_1
e_1 <- data.frame(bst_model_1$evaluation_log)
plot(e_1$iter, e_1$train_mlogloss, col = "blue")
lines(e_1$iter, e_1$test_mlogloss, col = "red")
min(e_1$test_mlogloss)
e_1[e_1$test_mlogloss == 0.18105, ] #Hence the 14th iteration was the lowest, though the lowest value has increased slightly

#In this case eta = 0.2 may not have been wise - stick with bst_model for now

#Feature Importance
imp <-  xgb.importance(colnames(dtrain), model = bst_model)
print(imp)
#Gain = improvement in accuracy by a feature to the branches it is on; most important column
xgb.plot.importance(imp)


#Prediction and confusin matrix - test data
#We will use the test data to create these because it is the most important data
p <- predict(bst_model, newdata = dtest)
head(p)
#The first patients probability is fairly high and therefore is most likely to have CVD? You add them as pairs and will get a probability of 1 - needs to be converted to use it in the confusion matrix
pred <- matrix(p, nrow = nc, ncol = length(p)/nc) %>%
  t() %>% #This will create a transpose for the matrix
  data.frame() %>%
  mutate(label = test_labels, max_prob = max.col(., "last")-1)
head(pred)
#"Max_prob" is what the model is predicting and "label" is what it is actually labelled as in the dataset
table(Prediction = pred$max_prob, Actual = pred$label)
#Seems like you haven't got the entire dataset as only 600 observations are shown, as opposed to 2000 observations

#More XGBoost Paramaters

#We can add more things onto the bst_model
bst_model_2 <- xgb.train(params = xgb_params,
                         data = dtrain,
                         nrounds = 100, #Maximum number of iterations
                         watchlist = watchlist,
                         eta = 0.2,
                         max.depth = 3,#Maximum tree depth and the default value is 6, can be increased or decreased to see if it improves the model; theoretically it can take values from 1 to infinity)
                         gamma = 0, #Default setting is actually 0, again can range from 1 to infinity; larger values > more conservative algorithm - you can increase gamma to reduce overfitting and to reduce the gap between the blue line and red line in the graph
                         subsample = 0.5, #Default value is 1, which is 100% and can take values from 0 to 1; lower values -> to prevent overfitting - if it is 0.5 for example the algorithm will randomly collect half of the data in stanses to grow the trees - to prevent overfitting
                         colsample_bytree = 1, #Default value is 1
                         missing = NA, #It guides missing values to left or right and then choose direction with higher gain with regard to the objective - useful when you have a lot of missing data
                         seed = 333) #Repeat results

#Now do the same as above but with the bst_model_2
bst_model_2
e_2 <- data.frame(bst_model_2$evaluation_log)
plot(e_2$iter, e_2$train_mlogloss, col = "blue")
lines(e_2$iter, e_2$test_mlogloss, col = "red")
min(e_2$test_mlogloss)
e[e_2$test_mlogloss == 0.173852, ] #19th iteration
imp_2 <-  xgb.importance(colnames(dtrain), model = bst_model_2)
print(imp_2)
xgb.plot.importance(imp_2)
p_2 <- predict(bst_model_2, newdata = dtest)
head(p_2)

pred_2 <- matrix(p_2, nrow = nc, ncol = length(p_2)/nc) %>%
  t() %>% #This will create a transpose for the matrix
  data.frame() %>%
  mutate(label = test_labels, max_prob = max.col(., "last")-1)
head(pred_2)

table(Prediction = pred_2$max_prob, Actual = pred_2$label)
