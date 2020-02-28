kNNImputeOptimization = function(data.in, log = 0, scaled = T,perParam = F, seed = 1234, plot = 0){
  #' This function will take in some data to impute (must be numeric)
  #' and a seed input. First, the complete cases of the data.in will 
  #' be selectd and put into the object data.complete. Next a random
  #' sample (of same size as number of rows there are in the complete set)
  #' of rows from the original data.in will be taken. A snapshot of the 
  #' missing data pattern (NApattern) of those rows from the original 
  #' data.in will be recorded. That same napatttern will be mapped into 
  #' the complete.cases. Now we have our original complete.cases and 
  #' a dataset with the same size but with a realistic data pattern that
  #' we have taken from the original (real missing data structure). 
  #' 
  #' Next we are going to perform knn imputation with several choices of k 
  #' record some sort of loss function between the real values that we know
  #' of and the imputed values with the different k. We are going to obtain
  #' a different loss function for each choice of k. 
  #' 
  #' @ scaled: (T or F) the code will scale your data and then impute it, the RMSE will be
  #' given for the scaled data as that is a more objetive way to visualise
  #' how well knn is doing (if features have different scales). 
  #' 
  #' @ perParam: (T or F) you can choose to have your MSE/RMSE per param as well as
  #'  a global one.
  #' 
  #'  The user can then run this algorithm as many times as he/she wants
  #'  to get a distribution of the error for each k. 
  
  require(impute)
  require(tidyverse)
  require(ggplot2)
  
  set.seed(seed) # this makes the sampling of different indices random
  # take complete data
  data.complete = data.in[complete.cases(data.in),]
  # within the range of rows from data.in grab as many rows as there are
  # rows in the complete data set
  rand.index = sample(x = 1:nrow(data.in), size = nrow(data.complete))
  # get the NA pattern in that random subset of the data.in
  NApattern = which(is.na(data.in[rand.index,]), arr.ind = T)
  rows = NApattern[,1]
  columns = NApattern[,2]
  
  # data to impute will be based in complete cases
  data.to.impute = data.complete
  # induce missingness in the data to impute with the real structure
  # of the original data -- feeding these indices does not seem to work
  # well
  # -- data.to.impute[rows, columns] = NA
  
  ## fill in NAs with nested for loop
  
  thePlaceToBe = environment()
  apply(NApattern, MARGIN=1,
        FUN=function(k){
          i=k[1];j=k[2]
          dat = get("data.to.impute", envir = thePlaceToBe)
          dat[i,j] = NA
          assign("data.to.impute",
                 dat, envir = thePlaceToBe)
        })
  
  #get mean and sds for scaling
  mean.cols = as.vector(colMeans(data.to.impute, na.rm = T))
  sd.cols = as.vector(apply(data.to.impute,2, function(col) sd(col, na.rm = T)))
  # scale the data to impute
  data.scaled = sweep(sweep(data.to.impute,
                            MARGIN = 2,mean.cols,'-'),
                      2,sd.cols,"/")
  # put in right format for imputation
  data.scaled = as.matrix(data.scaled)
  
  
  if (scaled == T){
    # -- if want scaled RMSE
    data.complete.scaled = sweep(sweep(data.complete,
                                       MARGIN = 2,mean.cols,'-'),
                                 2,sd.cols,"/")
    # put in right format for imputation
    data.complete.scaled = as.matrix(data.complete.scaled)
  }
  
  
  # BARE IN MIND THERE IS NO RANDOMNESS REGARDING impute.knn
  # it has a default random seed inside it
  # get list of imputed datasets for k = 1:20
  predictions.k = lapply(c(1:50),
                         function(x) 
                           #knn.impute(data.scaled, k = x, cat.var = NULL)
                           #knnImputation(data.scaled, k = x, scale = F)
                           impute.knn(data.scaled, k = x,
                                      rowmax = 1, colmax = 1)$data
  )
  
  # assert if there's any missing data in any of the imputed dataframes
  stopifnot(sum(sapply(1:length(predictions.k),
                       function(i) anyNA(predictions.k[[i]]))) == 0)
  # rescale data
  # -- if want scaled RMSE: predictions.k.rescaled = predictions.k 
  if (scaled == T){
    predictions.k.rescaled = predictions.k
  } else{
    predictions.k.rescaled = lapply(1:length(predictions.k),
                                    function(x)
                                      sweep(sweep(predictions.k[[x]],
                                                  MARGIN = 2,sd.cols,'*')
                                            ,2,mean.cols,"+"))
  }

  # get amount of missing data
  missing.dat = sum(is.na(data.scaled))
  
  ## MAYBE COMPUTE ERROR PER COLUMN AND GIVE BOXPLOT OF THAT
  ## try logging bio because of skewed variables
  # Calculate mean square error
  # apply for every dataset
  
  if (scaled == F){
    MSE = sapply(1:length(predictions.k),function(L){
      # mean of vector of squared differences between vector and original data
      mean(apply(NApattern, MARGIN=1,
                 FUN=function(k){
                   i=k[1];j=k[2]
                   # -- if want scaled RMSE
                   # (predictions.k.rescaled[[L]][i,j]-data.complete.scaled[i,j])**2
                   (predictions.k.rescaled[[L]][i,j]-data.complete[i,j])**2
                 }))
    })
  } else if (scaled == T){
    MSE = sapply(1:length(predictions.k),
                 function(L){
      # mean of vector of squared differences between vector and original data
      mean(apply(NApattern, MARGIN=1,
                 FUN=function(k){
                   i=k[1];j=k[2]
                   # -- if want scaled RMSE
                   (predictions.k.rescaled[[L]][i,j]-data.complete.scaled[i,j])**2
                 }))
    }
    #{mean(predictions.k.rescaled[[L]]-data.complete.scaled)**2} # this would take the means into account
    )
  } else{
    print("Argument for scaled must be TRUE or FALSE, please run your code again with a valid input")
    stopifnot(1 == 0)}
  
  if (perParam == T){ # get error per parameter for all k's by 
    # substracting predicted dataframe from original ones (all non 
    # inferreed values with cancel out) square and then take colmeans
    # then plot.
    MSEperParam = t(sapply(1:length(predictions.k),
                           function(L){
      # mean of vector of squared differences between vector and original data
      colSums((predictions.k.rescaled[[L]]-data.complete.scaled)**2)/colSums(is.na(data.scaled))
    }))
    
    # dataframe with 
    MSEperParam = data.frame(sqrt(MSEperParam))
    
    if (plot == 1){
      # could further implement to concatenate the 
      # per param results for several iterations,
      # though it's pretty infromative at the moment
      print("average RMSE per covariate (over 50 k values)")
      print(colMeans(MSEperParam))
      MSEperParam$k = 1:50
      print(MSEperParam %>% pivot_longer(-k) %>% 
              ggplot(aes(x = k, y = value))+
              facet_wrap(~name)+geom_point())
      # sort data.scaled in right order
      data.scaled = data.scaled[,sort(colnames(data.scaled))]
      
      print(MSEperParam %>% pivot_longer(-k) %>%
              group_by(name) %>% summarise(RMSE = mean(value)) %>% 
              mutate(RMSE = ifelse(is.na(RMSE), 0, RMSE)) %>%
              mutate(AmountNA = colMeans(is.na(data.scaled)))%>%
              ggplot(aes(x = reorder(name, RMSE), y = RMSE))+
              geom_col()+
              geom_text(aes(
                label = paste0(as.character(round(AmountNA,4)*100),
                               "%"), hjust = 1.2))+
              coord_flip())

      print(MSEperParam %>% pivot_longer(-k) %>%
              group_by(name) %>% summarise(RMSE = mean(value)) %>%
              mutate(RMSE = ifelse(is.na(RMSE), 0, RMSE)) %>% 
              mutate(AmountNA = colMeans(is.na(data.scaled)))%>%
              ggplot(aes(x = AmountNA, y = RMSE, label = name))+
              geom_point()+geom_text(aes(label = name), hjust = 0, vjust = 0)+
              ylim(0,2.5)+xlim(0,0.3))
    }

  }
  
  # calculate RMSE to have 
  RMSE = sqrt(MSE)
  RMSE
  
}