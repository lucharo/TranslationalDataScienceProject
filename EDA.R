# This script will load preprocessed datasets and plot the relevant EDA plots

##################################################################
##                     Input to this file:                      ##
##                      [bioProcessed.rds],                     ##
##                     [covProcessed.rds],                      ##
##                       [bioImputed.rds],
##                      [bioUnfiltered.rds]
##################################################################

#################################################################
##            Output from this file: many pdf plots            ##
#################################################################

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

library(parallel)
cores = detectCores()
##################################################################
##                  Load preprocessed datasets                  ##
##################################################################
bio.unfiltered = readRDS("data/preprocessed/bioUnfiltered.rds")
bio = readRDS("data/preprocessed/bioProcessed.rds")
cov = readRDS("data/preprocessed/covProcessed.rds")
bio.imp = readRDS("data/preprocessed/bioImputed.rds")

##################################################################
##                      Missing Data plots                      ##
##################################################################

# missing data plot with all 
bio.unfiltered.CVD = merge(bio.unfiltered,
                           cov[,"CVD_status"],
                           by="row.names",all.x=TRUE)

gg_miss_fct(x = bio.unfiltered.CVD, fct = CVD_status)+
  ggtitle("Missing data patterns by CVD status for all biomarkers")
ggsave("results/MissBioDataALLbyCVD.pdf")

bio.CVD = merge(bio,cov[,"CVD_status"],
                by="row.names",all.x=TRUE)

gg_miss_fct(x = bio.CVD, fct = CVD_status)+
  ggtitle("Missing data patterns by CVD status updated")
ggsave("results/MissBioDatabyCVD.pdf")

gg_miss_fct(x = cov, fct = CVD_status)

# There is no missing data in snp.original (checked full dataset on hpc)

# explanation of upset plot: the variables on the left are the set
# of variables with most missing data or whichever set of variables 
# you choose (nsets = n_var_miss(data), gives you all variables)
# with missing data, where as nsets = 10, gives you the top 10 vars
# with most missing data. The dots connected by lines represent
# the cases where NAs have been observed in the variables connected
# by the lines and on the top chart you have the number of ocurrences
# for each of those events. nintersects limits the amount of variable
# intersection you want to look at

gg_miss_upset(bio.unfiltered.CVD, 
                      nsets = 10,
                      nintersects = 10)
ggsave("results/upset_biofull_unfiltered.pdf")

upset = gg_miss_upset(bio.CVD, 
                      nsets = 10,
                      nintersects = 10)
ggsave("results/upset_biofull.pdf")

upset_cov = gg_miss_upset(cov)
ggsave("results/upset_covfull.pdf")

vis_miss(bio.unfiltered.CVD)+
  scale_y_continuous(position = 'right')+
  theme(axis.text.x = element_text(angle = 0))+
  scale_x_discrete(position = "bottom")+
  coord_flip()
ggsave("results/missBioDataPatterns_unfiltered.pdf")

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
b.pca = prcomp(bio.imp.CVD[,-ncol(bio.imp.CVD)], center = TRUE, scale. = TRUE)

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
