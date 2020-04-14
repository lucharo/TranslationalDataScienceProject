
# Load and Install relevant packages --------------------------------------

library(xgboost)
library(readr)
library(stringr)
library(caret)
library(car)
library(tidyverse)


# Merge Biomarkers and Covariates -----------------------------------------
bio = readRDS("data/preprocessed/bioImputed.rds")
cov = readRDS("data/preprocessed/covProcessed.rds")
bio.cov = merge(bio,cov,
                by="row.names",all.x=TRUE)
rownames(bio.cov) = bio.cov$Row.names
bio.cov = bio.cov[,-1]


# Setting up the environment ----------------------------------------------

#Set a random seed and shuffle data frame - doing this so that when I split our data into a tesing set and training set using the row numbers, I know that I'll get a random sample of data in both the testing and training set
set.seed(1234)
bio.cov <- bio.cov[sample(1:nrow(bio.cov)), ]
head(bio.cov)

# Remove information about the target variable from the training d --------

#First let's remove the columns that have information on our target variable
bio.cov_tar_rem <- bio.cov %>%
  select(-c("dc_cvd_st", "stop_cvd", "CVD_status", "age_CVD", "cvd_final_icd10", "primary_cause_death_ICD10", "cvd_death", "cvd_incident", "cvd_prevalent", "stop"))

#Let's create a new vector with labels - convert the CVD_status factor to an integer class starting at 0, as the first class should be 0; picky requirement
cov_label = as.integer(bio.cov$CVD_status)-1
nc_1 <- length(unique(cov_label))
nc_1

#Check out the first few lines
head(cov_label)

# Reduce the amount of redundant infromation ------------------------------

str(bio.cov_tar_rem)
summary(bio.cov_tar_rem)

#Select just the numeric columns
bio.cov_numeric <- bio.cov_tar_rem %>%
  select_if(is.numeric)
str(bio.cov_numeric)


# Convert dataframe into matrix -------------------------------------------

bio.cov_matrix <- data.matrix(bio.cov_numeric)
bio.cov_matrix
summary(bio.cov_matrix)

# Spliting Dataset into Training and Testing sets -------------------------

#Get the numb 70/30 training test split
numberofTrainingSamples_1 <- round(length(cov_label) * 0.7)

#Training Data
train_data_1 <- bio.cov_matrix[1:numberofTrainingSamples, ]
train_labels_1 <- cov_label[1:numberofTrainingSamples]

#Testing Data
test_data_1 <- bio.cov_matrix[-(1:numberofTrainingSamples), ]
test_labels_1 <- cov_label[-(1:numberofTrainingSamples)]


# Convert the cleaned dataframe into a DMATRIX ----------------------------

# put our testing & training data into two seperates Dmatrixs objects
dtrain_1 <- xgb.DMatrix(data = train_data_1, label= train_labels_1)
dtest_1 <- xgb.DMatrix(data = test_data_1, label= test_labels_1)


# Training our model ------------------------------------------------------

# train a model using our training data
model_1 <- xgboost(data = dtrain_1, # the data
                 nround = 2, # max number of boosting iterations
                 objective = "binary:logistic")

# generate predictions for our held-out testing data
pred_1 <- predict(model_1, dtest_1)

# get & print the classification error
err_1 <- mean(as.numeric(pred > 0.5) != test_labels_1)
print(paste("test-error=", err_1))


# Alternative Methods -----------------------------------------------------

xgb_params_1 <- list("objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = nc_1)

watchlist_1 <- list(train = dtrain_1, test = dtest_1)

#Extreme Gradeint Boosting Model
bst_model_b <- xgb.train(params = xgb_params_1,
                       data = dtrain_1,
                       nrounds = 100, #Maximum number of iterations
                       watchlist = watchlist_1 )


#Training and test error plot
bst_model_b
#Makes use of the evaluation_log results to make it into a plot
e_b <- data.frame(bst_model_b$evaluation_log)
plot(e_b$iter, e_b$train_mlogloss, col = "blue")

#Now let's add the test error to this plot
lines(e_b$iter, e_b$test_mlogloss, col = "red")


#To see the minimum value of the test-error
min(e_b$test_mlogloss)
#To know which iteration this minimum value belonged to
e_b[e_b$test_mlogloss == 0.151164, ] #Hence the 13th iteration was the lowest

#Now to optimise our model:
bst_model_b1 <- xgb.train(params = xgb_params_1,
                         data = dtrain_1,
                         nrounds = 100, #Maximum number of iterations
                         watchlist = watchlist_1,
                         eta = 0.2)
bst_model_b1
e_b1 <- data.frame(bst_model_b1$evaluation_log)
plot(e_b1$iter, e_b1$train_mlogloss, col = "blue")
lines(e_b1$iter, e_b1$test_mlogloss, col = "red")
min(e_b1$test_mlogloss)
e_b1[e_b1$test_mlogloss == 0.15154, ] #Hence the 17th iteration was the lowest

#Feature Importance
imp_b <-  xgb.importance(colnames(dtrain_1), model = bst_model_b1)
print(imp_b)
#Gain = improvement in accuracy by a feature to the branches it is on; most important column
xgb.plot.importance(imp_b)


#Prediction and confusin matrix - test data
#We will use the test data to create these because it is the most important data
p_b <- predict(bst_model_b1, newdata = dtest_1)
head(p_b)
#The first patients probability is fairly high and therefore is most likely to have CVD? You add them as pairs and will get a probability of 1 - needs to be converted to use it in the confusion matrix
pred_b <- matrix(p_b, nrow = nc_1, ncol = length(p_b)/nc_1) %>%
  t() %>% #This will create a transpose for the matrix
  data.frame() %>%
  mutate(label = test_labels_1, max_prob = max.col(., "last")-1)
head(pred_b)
#"Max_prob" is what the model is predicting and "label" is what it is actually labelled as in the dataset
table(Prediction = pred_b$max_prob, Actual = pred_b$label)
#Seems like you haven't got the entire dataset as only 600 observations are shown, as opposed to 2000 observations

#More XGBoost Paramaters

#We can add more things onto the bst_model
bst_model_2b <- xgb.train(params = xgb_params_1,
                         data = dtrain_1,
                         nrounds = 100, #Maximum number of iterations
                         watchlist = watchlist_1,
                         eta = 0.2,
                         max.depth = 3,#Maximum tree depth and the default value is 6, can be increased or decreased to see if it improves the model; theoretically it can take values from 1 to infinity)
                         gamma = 0, #Default setting is actually 0, again can range from 1 to infinity; larger values > more conservative algorithm - you can increase gamma to reduce overfitting and to reduce the gap between the blue line and red line in the graph
                         subsample = 0.5, #Default value is 1, which is 100% and can take values from 0 to 1; lower values -> to prevent overfitting - if it is 0.5 for example the algorithm will randomly collect half of the data in stanses to grow the trees - to prevent overfitting
                         colsample_bytree = 1, #Default value is 1
                         missing = NA, #It guides missing values to left or right and then choose direction with higher gain with regard to the objective - useful when you have a lot of missing data
                         seed = 333) #Repeat results
#######3 FROM OLD XGBCOV script






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





