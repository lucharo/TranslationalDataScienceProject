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

if (!require(ggbiplot)) install_github("vqv/ggbiplot")
library(ggbiplot)
if (!require(tidyverse)) install.packagaes("tidyverse")
library(tidyverse)
if (!require(naniar)) install.packages("naniar")
library(naniar)
if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if (!require(ggfortify)) install.packages("ggfortify")
library(ggfortify)
if (!require(dplyr)) install.packages("dplyr")
library(dplyr)
if (!require(stringr)) install.packages("stringr")
library(stringr)


library(parallel)
cores = detectCores()

#################################################################
##                     Setting environment                     ##
#################################################################
cluster = 0

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

ifelse(!dir.exists(file.path(save_plots, "EDA/")), dir.create(file.path(save_plots, "EDA/")), FALSE)
save_plots = paste0(save_plots,"EDA/")

#################################################################
##                       Saving function                       ##
#################################################################

save.results = function(figure, name, ggsv = T){
  if (ggsv){
    ggsave(paste0(save_plots,name,".pdf"))
  }
  saveRDS(figure, paste0(save_plots,name,".rds"))
}

##################################################################
##                  Load preprocessed datasets                  ##
##################################################################
bio.unfiltered = readRDS(paste0(data_folder,"bioUnfiltered.rds"))
bio = readRDS(paste0(data_folder,"bioProcessed.rds"))
cov = readRDS(paste0(data_folder,"covProcessed.rds"))
cov$CVD_status = as.factor(cov$CVD_status)
bio.imp = readRDS(paste0(data_folder,"bioImputedKNN.rds"))
snp = readRDS(paste0(data_folder,"snpProcessed.rds"))

# change dots for spaces in colnames bio for appropriate plotting
colnames(bio) = str_replace_all(colnames(bio), "\\."," ")



##################################################################
##                      Missing Data plots                      ##
##################################################################

# missing data plot with all (for biomarkers)
bio.unfiltered.CVD = merge(bio.unfiltered,
                           cov[,c("ID","CVD_status")],
                           by="ID")
# remove ID column as not needed for analysis
bio.unfiltered.CVD = bio.unfiltered.CVD[-1]

# gg_miss_fct(x = bio.unfiltered.CVD, fct = CVD_status)+
fig = bio.unfiltered.CVD %>%
  dplyr::group_by(CVD_status) %>%
  miss_var_summary() %>%
  ggplot(aes(CVD_status,
             variable,
             fill = pct_miss)) +
  geom_tile() +
  viridis::scale_fill_viridis(name = "% Miss") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle("Missing data patterns by CVD status for all biomarkers")

save.results(fig, "MissBioDataALLbyCVD")
# ggsave(paste0(save_plots,"MissBioDataALLbyCVD.pdf"))

bio.CVD = merge(bio,cov[,c("ID","CVD_status")],
                by="ID")
# remove ID col for plotting
bio.CVD = bio.CVD[-1]

# gg_miss_fct(x = bio.CVD, fct = CVD_status)+
fig = bio.CVD %>%
  dplyr::group_by(CVD_status) %>%
  miss_var_summary() %>%
  ggplot(aes(CVD_status,
             variable,
             fill = pct_miss)) +
  geom_tile() +
  viridis::scale_fill_viridis(name = "% Miss") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle("Missing data patterns by CVD status updated")
save.results(fig, "MissBioDatabyCVD")
# ggsave(paste0(save_plots,"MissBioDatabyCVD.pdf"))


# explanation of upset plot: the variables on the left are the set
# of variables with most missing data or whichever set of variables 
# you choose (nsets = n_var_miss(data), gives you all variables)
# with missing data, where as nsets = 10, gives you the top 10 vars
# with most missing data. The dots connected by lines represent
# the cases where NAs have been observed in the variables connected
# by the lines and on the top chart you have the number of ocurrences
# for each of those events. nintersects limits the amount of variable
# intersection you want to look at

# upset plots are not cool with ggsave method
pdf(file = paste0(save_plots,"upset_biofull_unfiltered.pdf"), onefile = F)
fig = gg_miss_upset(bio.unfiltered.CVD, 
                      nsets = 10,
                      nintersects = 10)
fig
dev.off()
save.results(fig, "upset_biofull_unfiltered", ggsv = F)


pdf(file = paste0(save_plots,"upset_biofull.pdf"), onefile = F)
fig = gg_miss_upset(bio.CVD, 
              nsets = 10,
              nintersects = 10)
fig
dev.off()
save.results(fig, "upset_biofull", ggsv = F)


miss1 = vis_miss(bio.unfiltered.CVD)+
  scale_y_continuous(position = 'right')+
  theme(axis.text.x = element_text(angle = 0))+
  scale_x_discrete(position = "bottom")+
  coord_flip()
ggsave(paste0(save_plots,"missBioDataPatterns_unfiltered.png")) # save as png because pdf creates artifact
saveRDS(miss1, paste0(save_plots,"missBioDataPatterns_unfiltered.rds"))

miss2 = vis_miss(bio.CVD)+
  scale_y_continuous(position = 'right')+
  theme(axis.text.x = element_text(angle = 0))+
  scale_x_discrete(position = "bottom")+
  coord_flip()
ggsave(paste0(save_plots,"missBioDataPatterns.png")) 
saveRDS(miss2, paste0(save_plots,"missBioDataPatterns.rds"))


# missing data plots for covariates

# gg_miss_fct(x = cov, fct = CVD_status) +
fig = cov %>%
  dplyr::group_by(CVD_status) %>%
  miss_var_summary() %>%
  ggplot(aes(CVD_status,
             variable,
             fill = pct_miss)) +
  geom_tile() +
  viridis::scale_fill_viridis(name = "% Miss") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle('Missing data patterns by CVD status for covariates')
# ggsave(paste0(save_plots,"MissCovDatabyCVD.pdf"))
save.results(fig, "MissCovDatabyCVD")

pdf(file = paste0(save_plots,"upset_covfull.pdf"))
fig = gg_miss_upset(cov)
fig
dev.off()
save.results(fig, "upset_covfull", ggsv = F)



# missing data plots for snps

snp.cvd = merge(snp, cov[,c("ID","CVD_status")], by="ID")
snp.cvd = snp.cvd[,-1]
# gg_miss_fct(x = snp.cvd, fct = CVD_status) +
snp.cvd %>%
  dplyr::group_by(CVD_status) %>%
  miss_var_summary() %>%
  ggplot(aes(CVD_status,
             variable,
             fill = pct_miss)) +
  geom_tile() +
  viridis::scale_fill_viridis(name = "% Miss") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle("Missing data patterns by CVD status for all snps")

snp = snp[,-1]
pdf(file = paste0(save_plots,"upset_snp.pdf"), onefile = F)
fig = gg_miss_upset(snp)
fig
dev.off()
save.results(fig, "upset_snp", ggsv = F)



##################################################################
##                   Biomarker distributions                    ##
##################################################################

bio.imp.CVD = merge(bio.imp,cov[,c("ID","CVD_status")],by="ID")
#rownames(bio.imp.CVD) = bio.imp.CVD$Row.names
bio.imp.CVD = bio.imp.CVD[,-1]

biomarker.labeller = function(original){
  str_replace_all(original, "\\."," ")
}


fig = bio.imp.CVD %>% pivot_longer(-CVD_status, 
                             names_to = "Biomarker",
                             values_to = "Amount") %>%
  ggplot(aes(x = Amount, color = CVD_status)) +
  geom_density(alpha = 0)+
  facet_wrap(Biomarker~.,
             scales = "free",
             labeller = labeller(Biomarker = biomarker.labeller))+
  #theme_minimal()+
  #scale_color_brewer(palette = "Set1")+
  ggtitle("Biomarker distribution for imputed data")

# ggsave(paste0(save_plots,"bio_dist_imp.pdf"))
save.results(fig, "bio_dist_imp")

fig = bio.imp.CVD %>% pivot_longer(-CVD_status, 
                             names_to = "Biomarker",
                             values_to = "Amount") %>%
  ggplot(aes(x = log10(Amount), color = CVD_status)) +
  geom_density(alpha = 0)+
  facet_wrap(Biomarker~.,
             scales = "free")+
  #theme_minimal()+
  #scale_color_brewer(palette = "Set1")+
  ggtitle("Biomarker distribution for imputed data(log scaled)")

# ggsave(paste0(save_plots,"bio_dist_impLOG.pdf"))
save.results(fig, "bio_dist_impLOG")

bio.CVD = bio.CVD[,-1]
fig = bio.CVD %>% pivot_longer(-CVD_status, 
                         names_to = "Biomarker",
                         values_to = "Amount") %>%
  ggplot(aes(x = Amount, color = CVD_status)) +
  geom_density(alpha = 0)+
  facet_wrap(Biomarker~., scales = "free")+
  #theme_minimal()+
  #scale_color_brewer(palette = "Set1")+
  #scale_colour_Publication()+
  ggtitle("Biomarker distribution for MCAR data")
# ggsave(paste0(save_plots,"bio_dist_MCAR.pdf"))
save.results(fig, "bio_dist_MCAR")

##################################################################
##                        PCA + PCA plot                        ##
##################################################################
# compute b's principal components,
# we set center and scale as TRUE, so that all features are scaled (~divided by sd) 
# and centred (~mean removed) before calculating PCAs as it should be.
b.pca = prcomp(bio.imp.CVD[,-ncol(bio.imp.CVD)],
               center = TRUE, scale. = TRUE)

# scree plot
scree = ggscreeplot(b.pca)+ylim(0,1)
# ggsave(paste0(save_plots,"Biomarkers_scree_plot.pdf"))
save.results(scree, "Biomarkers_scree_plot")

# ellipses with labels of CVD_status levels
ellipses = autoplot(b.pca, data = bio.imp.CVD,
                    colour = 'CVD_status',
                    alpha = 0.5,
                    loadings = T, loadings.colour = 'black',
                    loadings.label = T, loadings.label.colour = 'black')
# ggsave(paste0(save_plots,"PCA_plot_biomarkers.pdf"))
save.results(ellipses, "PCA_plot_biomarkers")


# print summary of pca
summary(b.pca)
saveRDS(b.pca, file = paste0(save_plots,"PCA_results_biomarkers.rds"))

# plotting original feature vectors and measurements projected onto top 2 PCAs
#ggbiplot(b.pca, groups = CVD_status)

