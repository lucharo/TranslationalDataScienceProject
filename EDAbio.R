#EDA

## NEXT STEPS:
#
# - apply to genetic data, which will be much more interesting
# - see if there is way to group thesemeasurements by disease outcome or something
# - investigate use of IDA and LDA
#
rm(list = ls())

###########################################################################
###########################################################################
###                                                                     ###
###                         PACKAGE DECLARATION                         ###
###                                                                     ###
###########################################################################
###########################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!require(devtools)) install.packages('devtools')
library(devtools)
if (!require(remotes)) install.packages('remotes')
library(remotes)
if (!require(ggbiplot)) install_github("vqv/ggbiplot")
library(ggbiplot)
if (!require(GGally)) install_github("GGally")
library(GGally)
if (!require(tidyverse)) install.packagaes("tidyverse")
library(tidyverse)
if (!require(naniar)) install.packages("naniar")
library(naniar)
if (!require(factoextra)) install.packages("factoextra")
library(factoextra)
if (!require(ggfortify)) install.packages("ggfortify")
library(ggfortify)
if (!require(stats)) install.packages("stats")
library(stats)
if (!require(mice)) install.packages('mice')
library(mice)


############################################################################
############################################################################
###                                                                      ###
###                             DATA LOADING                             ###
###                                                                      ###
############################################################################
############################################################################

#################################################################
##                      Original datasets                      ##
#################################################################
cov.original = readRDS("data/Covars_toy.rds")
bio.original= readRDS("data/Biomarkers_toy.rds")
bio.dict = readxl::read_xlsx("Biomarker_annotation.xlsx")
cov.dict = readxl::read_xlsx("Covariate_dictionary.xlsx")

##################################################################
##              Changing biomarkers codes by names              ##
##################################################################
#make nicely looking names (programmingly functional)
colnames(bio.dict) = make.names(colnames(bio.dict), unique=TRUE)

#get column numbers of columns with name containing pattern *(.)1(.)*
# use (.) to match the dot as opposed to using . as a wildcard
bio = bio.original[,!grepl("*(.)1(.)0", colnames(bio.original))]

# Match code with biomarker name to change column names of b
# get element 2 to 6 of all string in vector colnames(b)
# the match() function, match the substring from colnames to
# the UK.biobank.field in the biomarkers dictionary, 
# effectively ordering the colnames of b
# Alternative: order UK.bionbank.field entries and match them
#---- bio.dict = bio.dict %>% arrange(UK.Biobank.Field)

colnames(bio) = bio.dict$Biomarker.name[
  match(substring(colnames(bio),2,6),bio.dict$UK.Biobank.Field)] 

colnames(bio) = make.names(colnames(bio), unique=TRUE)
colnames(bio) = sub("\\.\\.",".", colnames(bio))

# safety-check for all vars being numeric
stopifnot(all(apply(bio, 2, is.numeric)))


##################################################################
##               Processing of covariates dataset               ##
##################################################################
## preprocessing c
# replace empty strings by NA
values_to_replace_w_na = c("")
cov = cov.original %>% replace_with_na_all(condition = ~.x %in%
                                             values_to_replace_w_na)
#remove anything to do wih cancer or external deaths
cov = cov[,!grepl("cancer|external", colnames(cov))]
#remove codes except for icd10
cov = cov[,!(colnames(cov) %in% c("cvd_final_icd9",
                                  "cvd_final_nc_illness_code",
                                  "cvd_final_opcs4",
                                  "cvd_final_ukb_oper_code",
                                  "other_cause_death"))]

# change numerical binary outcome variables to categorical
str(cov)
cov$CVD_status = as.factor(cov$CVD_status)
cov$vit_status = as.factor(cov$vit_status)
cov$dc_cvd_st = as.factor(cov$dc_cvd_st)
cov$cvd_death = as.factor(cov$cvd_death)



##################################################################
##                   Processing of biomarkers                   ##
##################################################################

# drop columns with more than 50% missing values i.e. Rheumatoid factor
# and oestradiol
bio = bio[,!colMeans(is.na(bio))>0.5]
saveRDS(bio, file = "data/bioProcessed.rds")

# impute biomarkers based on biomarkers only
t0 = Sys.time()
imp.model = mice(bio, m=5, maxit = 10,
                 seed = 500, printFlag = F)
print(Sys.time() - t0) # takes about 1 minute

# here we assign the imputed data to bio.imp
bio.imp = complete(imp.model,2)
saveRDS(bio.imp, file = "data/bioImputed.rds")

##################################################################
##                      Missing Data plots                      ##
##################################################################

bio.CVD = merge(bio,cov[,"CVD_status"],by="row.names",all.x=TRUE)

gg_miss_fct(x = bio.CVD, fct = CVD_status)+
  ggtitle("Missing data patterns by CVD status")
ggsave("results/MissBioDatabyCVD.pdf")

gg_miss_fct(x = cov, fct = CVD_status)

# explanation of upset plot: the variables on the left are the set
# of variables with most missing data or whichever set of variables 
# you choose (nsets = n_var_miss(data), gives you all variables)
# with missing data, where as nsets = 10, gives you the top 10 vars
# with most missing data. The dots connected by lines represent
# the cases where NAs have been observed in the variables connected
# by the lines and on the top chart you have the number of ocurrences
# for each of those events. nintersects limits the amount of variable
# intersection you want to look at

upset = gg_miss_upset(bio.CVD, 
                      nsets = 10,
                      nintersects = 10)


vis_miss(bio.CVD)+
  scale_y_continuous(position = 'right')+
  theme(axis.text.x = element_text(angle = 0))+
  scale_x_discrete(position = "bottom")+
  coord_flip()
ggsave("results/missBioDataPatterns.pdf")

##################################################################
##                   Biomarker distributions                    ##
##################################################################

bio.imp.CVD = merge(bio.imp,cov[,"CVD_status"],by="row.names",all.x=TRUE)
rownames(bio.imp.CVD) = bio.imp.CVD$Row.names
bio.imp.CVD = bio.imp.CVD[,-1]
bio.imp.CVD %>% pivot_longer(-CVD_status, 
                           names_to = "Biomarker",
                           values_to = "Amount") %>%
  ggplot(aes(x = Amount, color = CVD_status)) +
  geom_density(alpha = 0)+
  facet_wrap(Biomarker~., scales = "free")+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")+
  ggtitle("Biomarker distribution for imputed data")
ggsave("results/bio_dist_imp.pdf")

rownames(bio.CVD) = bio.CVD$Row.names
bio.CVD = bio.CVD[,-1]
bio.CVD %>% pivot_longer(-CVD_status, 
                         names_to = "Biomarker",
                         values_to = "Amount") %>%
  ggplot(aes(x = Amount, color = CVD_status)) +
  geom_density(alpha = 0)+
  facet_wrap(Biomarker~., scales = "free")+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")+
  ggtitle("Biomarker distribution for MCAR data")
ggsave("results/bio_dist_MCAR.pdf")

##################################################################
##                        PCA + PCA plot                        ##
##################################################################
# compute b's principal components,
# we set center and scale as TRUE, so that all features are scaled (~divided by sd) 
# and centred (~mean removed) before calculating PCAs as it should be.
b.pca = prcomp(bio.imp.CVD[,-28], center = TRUE, scale. = TRUE)

# scree plot
scree = ggscreeplot(b.pca)+ylim(0,1)
ggsave("results/Biomarkers_scree_plot.pdf")

# ellipses with labels of CVD_status levels
ellipses = autoplot(b.pca, data = bio.imp.CVD,
                    colour = 'CVD_status',
                    alpha = 0.5,
                    loadings = T, loadings.colour = 'black',
                    loadings.label = T, loadings.label.colour = 'black')+
  scale_color_brewer(palette = 'Set1')
ggsave("results/PCA_plot_biomarkers.pdf")

# print summary of pca
summary(b.pca)
save(b.pca, file = "results/PCA_results_biomarkers.rds")

# plotting original feature vectors and measurements projected onto top 2 PCAs
#ggbiplot(b.pca, groups = CVD_status)






##################################################################
##                         ggpairs plot                         ##
##################################################################

#get "interesting" biomarkers name from the non-complete dataset to plot
# ggpairs of those

# interesting_biomarkers = colnames(readRDS("data/Biomarkers_toy.rds"))
# interesting_biomarkers = sub(".0.0", "", interesting_biomarkers)
# interesting_biomarkers = gsub("_", ".", interesting_biomarkers)

# this would be good but the biomarkers do not even match

# pairs plot of biomarkers separated by vital status

#ggpairs(bio.imp,
        #aes(color=CVD_status), progress = FALSE)




##################################################################
##                         First models                         ##
##################################################################
## Running basic models
#Covariates
# glm_cov <- glm(CVD_status ~ age_cl + BS2_all, data=cov, family=binomial)
# summary(glm_cov)
# round(cbind("odds" = exp(coef(glm_cov)), exp(confint(glm_cov))), 3)
# 
# if (!require(jtools)) install.packages('jtools')
# library(jtools)
# glm_cov_plot <- effect_plot(glm_cov, pred = BS2_all, interval = TRUE, int.width = 0.95) 
# glm_cov_plot
# 
# #Biomarkers
# # all_cols <- colnames(bio.joint[1:30]) ---- this line's not needed
# glm_bio <- glm(CVD_status ~., data=bio.joint, family=binomial)
# summary(glm_bio)
# round(cbind("odds" = exp(coef(glm_bio)), exp(confint(glm_bio))), 3)















