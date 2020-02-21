
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





