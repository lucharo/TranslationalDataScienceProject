## PLS/Elastic Net attempt on TDS data


library(glmnet)
library(lme4)
# Aim of this script is to replicate the work from practical 3 and 4 on the
# TDS dataset

# Maybe also do univariate analysis using mixed models effect (maybe for
# duplicate measurements)

# load biomarkers dataset, ideally one of the already preprocessed ones
# so you dont have to run preprocessing again


###########################################################################
###########################################################################
###                                                                     ###
###                UNIVARIATE ANALYSIS OF THE BIOMARKERS                ###
###                                                                     ###
###########################################################################
###########################################################################

bio.covs = readRDS("data/bio_and_covs.rds")

confounders = c("age_CVD","BS2_all", "qual2","smok_ever_2",
                "physical_activity_2", "alcohol_2", "BMI_5cl_2",
                "no_cmrbt_cl2", "no_medicines")

bio.names = colnames(select(biomarkersCVD, -CVD_status))

DoYouMatter = function(molecule, data = bio.covs){
  # I guess I would ideally not have to precise a data argument and allow
  # the data to be the global proteins.covars object, though lmer throws an
  # error
  
  formula0 = as.formula(paste(molecule,
                              paste(confounders, collapse = "+"), sep = "~"))
  m0 = lm(formula0, data = data, na.action = na.omit)
  
  formula1 = as.formula(paste(molecule,
                              paste(c(confounders, "CVD_status"),
                                    collapse = "+"),
                              sep = "~"))
  m1 = lm(formula1, data = data, na.action = na.omit)
  
  anova(m0,m1)$`Pr(>Chisq)`[2]
}

pvals = unlist(lapply(bio.names, DoYouMatter))
FatANOVA = data.frame("Biomarkers" = bio.names, "p-value" = pvals)

