library(tidyverse)
library(ggpubr)

cluster = 0
platform = Sys.info()['sysname']
if (platform=='Linux'){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

ifelse(!dir.exists(file.path(save_plots, "BHS/")), dir.create(file.path(save_plots, "BHS/")), FALSE)
save_plots = paste0(save_plots,"BHS/")

#### ADD COMPARISONS BY SUBGROUP and add some more comments

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))
bio.imp = readRDS(paste0(data_folder,"bioImputedKNN.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))

# bio.dict = readxl::read_xlsx("../Biomarker_annotation.xlsx")
#
# # make pretty names
# bio.dict$`Biomarker name` = make.names(bio.dict$`Biomarker name`)
# bio.dict$`Biomarker name` = sub("\\.\\.",".",
#                                 bio.dict$`Biomarker name`)

# for simplicity we merge cov and bio here
bio.imp = merge(bio.imp, cov[,c("ID","age_cl","gender")], by = "ID")
ids = bio.imp$ID # wanna keep an explicit copy od the IDs
rownames(bio.imp) = bio.imp[,1]
bio.imp = bio.imp[,-1]

bio.imp$age_cl = as.factor(bio.imp$age_cl)
bio.imp$gender = as.factor(bio.imp$gender)
#
BHSCalculator = function(reference, stratified = F, bySystems = T){
  # reference can take values:
  # -- "More_is_bad_Mantej
  # or
  # -- "More_is_bad_paper

  MoreIsBad = paste0("MoreIsBad.",reference)

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
                                    AgeClass = character(0),
                                    Gender = character(0))

    for (gender in unique(bio.imp$gender)){
      for (age.class in unique(bio.imp$age_cl)){

        # getting the entries from the bio dataset where the cov$age_cl values are
        # equal to age.class
        bio.per.class = bio.imp[bio.imp$age_cl==age.class & bio.imp$gender == gender,
                                -c(ncol(bio.imp)-1, ncol(bio.imp))]
        # -ncol(bio.imp) everywhere to avoid age_cl column

        result = data.frame(Biomarker = colnames(bio.per.class),
                            data.frame(matrix(
                              unlist(
                                lapply(colnames(bio.per.class),
                                       function(x) quantile_check(x, reference = MoreIsBad, bio.per.class))),
                              nrow=length(colnames(bio.per.class)),
                              byrow=T), stringsAsFactors = F),
                            stringsAsFactors = F)

        colnames(result)[2:3] = c("Quantile.value", "Quartile")
        result = cbind(result, AgeClass = age.class, Gender = gender)

        relevant_quantiles = rbind(relevant_quantiles, result)
      }
    }

    # take only complete cases
    relevant_quantiles = relevant_quantiles[complete.cases(relevant_quantiles),]

  } else {
    relevant_quantiles = data.frame(Biomarker = colnames(bio.imp[,-c(ncol(bio.imp)-1, ncol(bio.imp))]),
                                    data.frame(matrix(
                                      unlist(
                                        lapply(colnames(bio.imp[,-c(ncol(bio.imp)-1, ncol(bio.imp))]),
                                               function(x) quantile_check(x, reference = MoreIsBad, bio.imp))),
                                      nrow=length(colnames(bio.imp[,-c(ncol(bio.imp)-1, ncol(bio.imp))])),
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

    # initialising empty dataframe with dims of however many biomarkers are being considered
    # could add two extra columns for each strata
    bio.score = bio.imp[F,unique(relevant_quantiles$Biomarker)]

    for (gender in unique(bio.imp$gender)){
      for (age in unique(bio.imp$age_cl)){
        # select entries of biomarker for which age is equal to iterations
        bio.sub  = bio.imp[bio.imp$age_cl==age & bio.imp$gender == gender,
                           -c(ncol(bio.imp)-1, ncol(bio.imp))]
        # select quantiles for which age is equal to age iterations
        relevant_quantiles.age.gender = relevant_quantiles[relevant_quantiles$AgeClass==age &
                                                             relevant_quantiles$Gender == gender,]
        attach(relevant_quantiles.age.gender)
        res = lapply(Biomarker, # for each biomarker from relevant_quantiles.age.gender
                     function(x){ # if the biomarker is bad in excess (i.e. boundary is 3rd quartile)
                       # return the element-wise comparison to the relevant quartile value
                       ifelse(Quartile[Biomarker == x] == "3rd",
                              return(bio.sub[,x] >
                                       Quantile.value[Biomarker == x]),
                              return(bio.sub[,x] <
                                       Quantile.value[Biomarker == x]))
                     })
        res = data.frame(matrix(unlist(res),
                                ncol = length(Biomarker), byrow = F), # put in right format
                         stringsAsFactors = F)
        rownames(res) = rownames(bio.sub)
        colnames(res) = Biomarker
        #res$AgeCl = age
        bio.score = rbind(bio.score, res)

        # stopifnot(all(bio.imp$age_cl[bio.imp$age_cl==age]==res$AgeCl))
        # stopifnot(all(rownames(bio.imp$age_cl[bio.imp$age_cl==age])==rownames(res$AgeCl)))

        detach(relevant_quantiles.age.gender)
      }
    }


    # I actually don't need the age class column just need them to keep the original index
    #bio.score = bio.score[,-c(ncol(bio.score)-1, ncol(bio.score))]

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
  numbersofBio = ncol(bio.score)
  ##################################################################
  ##    3rd and finally get the mean score for each individual    ##
  ##              and sort in same order as bio.imp               ##
  ##################################################################


  if (bySystems){
    systems.ref = paste0("System.",reference)
    bio.score$ID = rownames(bio.score)
    bio.score = bio.score %>% tidyr::gather(key = "Biomarker", value = "Amount", -ID) %>%
      merge(y = bio.dict[,c(systems.ref, "Biomarker name")], all.x = T, by.x = "Biomarker", by.y = "Biomarker name")
    colnames(bio.score)[ncol(bio.score)] = "System"
    # get score for each individual by ID and individual system
    bio.score = bio.score %>% dplyr::group_by(ID, System) %>% dplyr::summarise(BHSsystem = mean(Amount))
    # the line below would give you dataframe with score by system
    # bio.score %>% spread(System, BHS)

    bio.score = bio.score %>% dplyr::group_by(ID) %>% dplyr::summarise(total_score = mean(BHSsystem))
    bio.score = as.data.frame(bio.score) # Setting row names on a tibble is deprecated.
    rownames(bio.score) = bio.score$ID
  }else {
    # get average "grade" over all the biomarkers
    bio.score$total_score = rowMeans(bio.score)
  }

  # sort bio.score in by the row.name as bio.imp
  bio.score = bio.score[rownames(bio.imp),]
  stopifnot(all(rownames(bio.imp)==rownames(bio.score)))
  print(paste(numbersofBio," biomarkers where assessed for this BHS calculation"))
  bio.score$total_score
}

scores_paper = BHSCalculator("Paper",T)
# at the end of the BHS calculator I have a function that ensures
# bio.imp and bio.score are in the same order
ScoresPaper = as.data.frame(
  cbind(ID = ids, BHS = scores_paper), stringsAsFactors = F)
saveRDS(ScoresPaper, paste0(save_plots, "ScoresPaper.rds"))

scores_Mantej = BHSCalculator("Mantej",T)
ScoresMantej = as.data.frame(
  cbind(ID = ids, BHS = scores_Mantej), stringsAsFactors = F)
saveRDS(ScoresMantej, paste0(save_plots, "ScoresMantej.rds"))

scores_Barbara = BHSCalculator("Barbara",T)
ScoresBarbara = as.data.frame(
  cbind(ID = ids, BHS = scores_Barbara), stringsAsFactors = F)
saveRDS(ScoresBarbara, paste0(save_plots, "ScoresBarbara.rds"))

##################################################################
##                             Plot                             ##
##################################################################

my_comparisons <- list(c("BS2", "Barbara"), c("Mantej", "Paper"),c("BS2", "Mantej"),
                       c("Barbara", "Mantej"), c("BS2", "Paper"),c("Barbara", "Paper"))

fig = data.frame(rbind(
  cbind(BHS = scores_paper, Reference = "Paper"),
  cbind(BHS = scores_Mantej, Reference = "Mantej"),
  cbind(BHS = scores_Barbara, Reference = "Barbara"),
  cbind(BHS = cov$BS2_all[cov$ID %in% ScoresBarbara$ID], Reference = "BS2")),
      stringsAsFactors = F) %>%
  ggplot(aes(x = Reference, y = as.numeric(BHS), fill = Reference))+
  geom_boxplot()+
  ylab("BHS")+
  scale_fill_brewer(palette = "Set1")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test",
                     paired = F, label.y = c(1,1,1,1.05, 1.1, 1.15), label = "p.signif")+
  stat_summary(geom = "point", shape = 23)+
  theme_minimal()+ggtitle("Distribution of the BHS Scores")
fig

ggsave(paste0(save_plots, "MantejPaperBarbara.pdf"))
saveRDS(fig, paste0(save_plots, "MantejPaperBarbara.rds"))


##################################################################
##                            t-test                            ##
##################################################################
# print(t.test(as.numeric(BHS) ~ Reference, data = data.frame(
#   rbind(cbind(BHS = scores_paper, Reference = "Paper"),
#         cbind(BHS = scores_Mantej, Reference = "Mantej"),
#         cbind(BHS = scores_Barbara, Reference = "Barbara")),
#   stringsAsFactors = F), paired = T))

jointdata = data.frame(
  rbind(
    cbind(BHS = scores_paper, Reference = "Paper"),
    cbind(BHS = scores_Mantej, Reference = "Mantej"),
    cbind(BHS = scores_Barbara, Reference = "Barbara"),
    cbind(BHS = cov$BS2_all[cov$ID %in% ScoresBarbara$ID], Reference = "BS2")
  )
)
jointdata$BHS = as.numeric(jointdata$BHS)

print(compare_means(BHS ~ Reference,
              method = "t.test",
              paired = T,
              data = jointdata))


## BasicModel
Scores.CVD.Mantej = merge(ScoresMantej,
                          cov[, c("CVD_status","ID")], by = "ID")
Scores.CVD.Paper = merge(ScoresPaper,
                         cov[, c("CVD_status","ID")], by = "ID")
Scores.CVD.Barbara = merge(ScoresBarbara,
                         cov[, c("CVD_status","ID")], by = "ID")

Scores.CVD.Mantej$CVD_status = as.factor(Scores.CVD.Mantej$CVD_status)
Scores.CVD.Paper$CVD_status = as.factor(Scores.CVD.Paper$CVD_status)
Scores.CVD.Barbara$CVD_status = as.factor(Scores.CVD.Barbara$CVD_status)

mod1 = glm(CVD_status~BHS, data = Scores.CVD.Mantej, family = "binomial")
print("Mantej score:")
print(mod1)
print(exp(coef(mod1)))
print(exp(confint(mod1)))
print(summary(mod1)$coefficients)

print("##################################################################")

mod2 = glm(CVD_status~BHS, data = Scores.CVD.Paper, family = "binomial")
print("Paper score:")
print(mod2)
print(exp(coef(mod2)))
print(exp(confint(mod2)))
print(summary(mod2)$coefficients)

print("##################################################################")

mod3 = glm(CVD_status~BHS, data = Scores.CVD.Barbara, family = "binomial")
print("Barbara score:")
print(mod3)
print(exp(coef(mod3)))
print(exp(confint(mod3)))
print(summary(mod3)$coefficients)

print("##################################################################")

mod4 = glm(as.factor(CVD_status)~BS2_all, data = cov, family = "binomial")
print("BS2_all score:")
print(mod4)
print(exp(coef(mod4)))
print(exp(confint(mod4)))
print(summary(mod4)$coefficients)

