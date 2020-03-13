library(tidyverse)
library(ggpubr)

cluster = 0

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

#### ADD COMPARISONS BY SUBGROUP and add some more comments

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))
bio.imp = readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))

bio.dict = readxl::read_xlsx("../Biomarker_annotation.xlsx")

# make pretty names
bio.dict$`Biomarker name` = make.names(bio.dict$`Biomarker name`)
bio.dict$`Biomarker name` = sub("\\.\\.",".",
                                bio.dict$`Biomarker name`)

# for simplicity we merge cov and bio here
bio.imp = merge(bio.imp, cov[,c("age_cl")], by = "row.names")
rownames(bio.imp) = bio.imp[,1]
bio.imp = bio.imp[,-1]

BHSCalculator = function(reference, stratified = F){
  # reference can take values: 
  # -- "More_is_bad_Mantej
  # or
  # -- "More_is_bad_paper
  
  quantile_check = function(column, reference, dataset){
    # the "reference" column in this function is the one containing the info on whether or not
    # a given molecule is harmful in excess or in shortage, there are two references:
    # -- One from the biomarkers available in the paper
    # -- Another one created from expert knowledge (by Mantej) including more biomarkers than 
    # in the original model
    
    # if the value of this reference column for a given biomarker is 1, it means that biomarker
    # is harmful in excess, in which case we select the 3rd quartile (75%) which will be our benchmark
    # MORE than this kind of benchmark ==> +1 BHS
    
    if (bio.dict[bio.dict$`Biomarker name` == column,
                 reference] == 1){
        output = quantile(dataset[,column])[4]
        quartile = "3rd"
    }
    # in the case the value of the reference column for a given biomarker is 0, it means the biomarker 
    # is harmful in shortage, in which case we select the 1st quartile (25%) which will be our bencmark
    # LESS than this kind of benchmark ==> +1 BHS
    else if (bio.dict[bio.dict$`Biomarker name` == column,
                      reference] == 0){
        output = quantile(dataset[,column])[2]
        quartile = "1st"
      }
    # Store info on biomarker with missing reference value, still will be removed later
    # In the references provided by Mantej, we wrote -1 for those biomarkers we were unsure if they 
    # were beneficial or harmful (or irrelevant)
    #
    else{
      output = NA
      quartile = NA
    }
    list(output, quartile)
  }
  
  
  ##################################################################
  ##        First make a 'dictionary' storing the relevant        ##
  ##         quantiles for each biomarker and each strata         ##
  ##################################################################
  # store relevant quartile info in a dataframe
  if (stratified){
    relevant_quantiles = data.frame(Biomarker = character(0),
                                    Quantile.value = numeric(0),
                                    Quartile = character(0),
                                    AgeClass = numeric(0))
    for (age.class in 1:3){
      
      # getting the entries from the bio dataset where the cov$age_cl values are 
      # equal to age.class
      bio.per.class = bio.imp[bio.imp$age_cl==age.class,-ncol(bio.imp)]
                                        # -ncol(bio.imp) everywhere to avoid age_cl column
      
      result = data.frame(Biomarker = colnames(bio.per.class),
                          data.frame(matrix(
                            unlist(
                              lapply(colnames(bio.per.class),
                                     function(x) quantile_check(x, reference = reference, bio.per.class))),
                            nrow=length(colnames(bio.per.class)),
                            byrow=T), stringsAsFactors = F), 
                          stringsAsFactors = F)
      
      colnames(result)[2:3] = c("Quantile.value", "Quartile")
      result = cbind(result, AgeClass = age.class)
      
      relevant_quantiles = rbind(relevant_quantiles, result)
    }
    # take only complete cases
    relevant_quantiles = relevant_quantiles[complete.cases(relevant_quantiles),]
    
  } else {
    relevant_quantiles = data.frame(Biomarker = colnames(bio.imp[,-ncol(bio.imp)]),
                                    data.frame(matrix(
                                      unlist(
                                        lapply(colnames(bio.imp[,-ncol(bio.imp)]),
                                               function(x) quantile_check(x, reference = reference, bio.imp))),
                                      nrow=length(colnames(bio.imp[,-ncol(bio.imp)])),
                                      byrow=T), stringsAsFactors = F), 
                                    stringsAsFactors = F)
    
    
    colnames(relevant_quantiles)[2:3] = c("Quantile.value", "Quartile")
    
    # take only complete cases
    relevant_quantiles = relevant_quantiles[complete.cases(relevant_quantiles),]
  }

  
  #######################################################################
  ##  Second go through every biomarker (and each strata if relevant)  ##
  ##                and get the score for each biomarker               ##
  #######################################################################
  # bio.scores stores for each individual whether its values are over/under the relevant quantiles
  # compare individual values for each biomark to their "benchmark"
  if (stratified){

    # initialising empty dataframe with dims of bio.imp
    bio.score = bio.imp[F,1:(length(unique(relevant_quantiles$Biomarker))+1)]
    
    for (age in 1:3){
      # select entries of biomarker for which age is equal to iterations
      bio.sub  = bio.imp[bio.imp$age_cl==age,-ncol(bio.imp)]
      # select quantiles for which age is equal to age iterations
      relevant_quantiles.age = relevant_quantiles[relevant_quantiles$AgeClass==age,]
      attach(relevant_quantiles.age)
      res = lapply(Biomarker,
                         function(x){
                           ifelse(Quartile[Biomarker == x] == "3rd",
                                  return(bio.sub[,x] > 
                                           Quantile.value[Biomarker == x]),
                                  return(bio.sub[,x] < 
                                           Quantile.value[Biomarker == x]))
                         })
      res = data.frame(matrix(unlist(res),
                                    ncol = length(Biomarker), byrow = F),
                             stringsAsFactors = F)
      rownames(res) = rownames(bio.sub)
      colnames(res) = Biomarker
      res$AgeCl = age
      bio.score = rbind(bio.score, res)
      
      stopifnot(all(bio.imp$age_cl[bio.imp$age_cl==age]==res$AgeCl))
      stopifnot(all(rownames(bio.imp$age_cl[bio.imp$age_cl==age])==rownames(res$AgeCl)))
      
      detach(relevant_quantiles.age)
    }
  
    # I actually don't need the age class column just need them to keep the original index
    bio.score = bio.score[,-ncol(bio.score)]
 
  }else{
    attach(relevant_quantiles)
    # bio.score is a list containing for each biomarker (first subset of list)
    # for each individual whether the value of that biomarker is in the healthy/unhealthy range(FALSE, TRUE)
    # respectively
    bio.score = lapply(Biomarker,
                       function(x){
                         ifelse(Quartile[Biomarker == x] == "3rd",
                                return(bio.imp[,x] > Quantile.value[Biomarker == x]),
                                return(bio.imp[,x] < Quantile.value[Biomarker == x]))
                       })
    
    # bio.score is then turned into a nice data frame where each biomarker is a column
    # and each row is an individual 
    # note: unlist returns a long vector which results from unlist the bio.score list into
    # a vector with a length equal to the number of elements in the original list.
    # byrow in this case must be equal to FALSE
    bio.score = data.frame(matrix(unlist(bio.score),
                                  ncol = length(Biomarker), byrow = F),
                           stringsAsFactors = F)
    rownames(bio.score) = rownames(bio.imp)
    detach(relevant_quantiles)
  }

  # This is going on for every row (in code above)
  # bio.imp$Testosterone < relevant_quantiles$Quantile.value[relevant_quantiles$Biomarker == "Testosterone"]
  
  ##################################################################
  ##    3rd and finally get the mean score for each individual    ##
  ##              and sort in same order as bio.imp               ##
  ##################################################################
  # get average "grade" over all the biomarkers
  bio.score$total_score = rowMeans(bio.score)
  
  # sort bio.score in the same order as bio.imp
  bio.score = bio.score[rownames(bio.imp),]
  stopifnot(all(rownames(bio.imp)==rownames(bio.score)))
  print(paste(length(colnames(bio.score))," biomarkers where assessed for this BHS calculation"))
  bio.score$total_score
}

scores_paper = BHSCalculator("More_is_bad_paper",T)

scores_Mantej = BHSCalculator("More_is_bad_Mantej",T)




##################################################################
##                             Plot                             ##
##################################################################
data.frame(rbind(cbind(BHS = scores_paper, Reference = "Paper"),
      cbind(BHS = scores_Mantej, Reference = "Mantej")),
      stringsAsFactors = F) %>% 
  ggplot(aes(x = Reference, y = as.numeric(BHS), fill = Reference))+
  geom_boxplot()+
  ylab("BHS")+
  scale_fill_brewer(palette = "Set1")+
  stat_compare_means(method = "t.test", paired = T,
                     label.x = 1.5, label.y = 0.9)+
  stat_summary(geom = "point", shape = 23)+
  theme_minimal()


##################################################################
##                            t-test                            ##
##################################################################
t.test(as.numeric(BHS) ~ Reference, data = data.frame(
  rbind(cbind(BHS = scores_paper, Reference = "Paper"),
        cbind(BHS = scores_Mantej, Reference = "Mantej")),
  stringsAsFactors = F), paired = T)
