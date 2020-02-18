library(tidyverse)
library(ggpubr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

bio = readRDS("data/preprocessed/bioProcessed.rds")
bio.imp = readRDS("data/preprocessed/bioImputed.rds")
cov = readRDS("data/preprocessed/covProcessed.rds")

bio.dict = readxl::read_xlsx("Biomarker_annotation.xlsx")

bio.dict$`Biomarker name` = make.names(bio.dict$`Biomarker name`)
bio.dict$`Biomarker name` = sub("\\.\\.",".",
                                bio.dict$`Biomarker name`)



BHSCalculator = function(reference){
  # reference can take values: 
  # -- "More_is_bad_Mantej
  # or
  # -- "More_is_bad_paper
  
  quantile_check = function(column, reference){
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
      output = quantile(bio.imp[,column])[4]
      quartile = "3rd"
    }
    # in the case the value of the reference column for a given biomarker is 0, it means the biomarker 
    # is harmful in shortage, in which case we select the 1st quartile (25%) which will be our bencmark
    # LESS than this kind of benchmark ==> +1 BHS
    else if (bio.dict[bio.dict$`Biomarker name` == column,
                      reference] == 0){
      output = quantile(bio.imp[,column])[2]
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
  
  # store relevant quartile info in a dataframe
  relevant_quantiles = data.frame(Biomarker = colnames(bio.imp),
                                  data.frame(matrix(
                                    unlist(
                                      lapply(colnames(bio.imp),
                                             function(x) quantile_check(x, reference = reference))),
                                    nrow=length(colnames(bio.imp)),
                                    byrow=T), stringsAsFactors = F), 
                                  stringsAsFactors = F)
  
  colnames(relevant_quantiles)[2:3] = c("Quantile.value", "Quartile")
  
  # take only complete cases
  relevant_quantiles = relevant_quantiles[complete.cases(relevant_quantiles),]
  
  # bio.scores stores for each individual whether its values are over/under the relevant quantiles
  # compare individual values for each biomark to their "benchmark"
  attach(relevant_quantiles)
  bio.score = lapply(Biomarker,
                     function(x){
                       ifelse(Quartile[Biomarker == x] == "3rd",
                              return(bio.imp[,x] > Quantile.value[Biomarker == x]),
                              return(bio.imp[,x] < Quantile.value[Biomarker == x]))
                     })
  bio.score = data.frame(matrix(unlist(bio.score),
                                ncol = length(Biomarker), byrow = T),
                         stringsAsFactors = F)
  detach(relevant_quantiles)
  
  
  # This is going on for every row
  # bio.imp$Testosterone < relevant_quantiles$Quantile.value[relevant_quantiles$Biomarker == "Testosterone"]
  
  bio.score$total_score = rowSums(bio.score)/length(colnames(bio.score))
  print(paste(length(colnames(bio.score))," biomarkers where assessed for this BHS calculation"))
  bio.score$total_score
}

scores_paper = BHSCalculator("More_is_bad_paper")

scores_Mantej = BHSCalculator("More_is_bad_Mantej")

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


t.test(as.numeric(BHS) ~ Reference, data = data.frame(
  rbind(cbind(BHS = scores_paper, Reference = "Paper"),
        cbind(BHS = scores_Mantej, Reference = "Mantej")),
  stringsAsFactors = F), paired = T)
