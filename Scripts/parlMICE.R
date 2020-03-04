
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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
if (!require(DMwR)) install.packages('DMwR')
library(DMwR)
if (!require(impute)) BiocManager::install("impute")
library(impute)
library(parallel)

cluster = 0
save_folder = data_folder = "../data/preprocessed/"
if (cluster == 1){
  save_folder = data_folder = "../FULLDATA/preprocessed/"
}

bio = readRDS(paste0(data_folder,"bioProcessed.rds"))

# Impute with parallelisation
t0 = Sys.time()
# parlmice handles how many cores to set by itself
imp.model = parlmice(bio, m=5, seed = NA, printFlag = F,
                       cl.type = "FORK")
print(Sys.time() - t0) # takes about 1 minute
# 
# # here we assign the imputed data to bio.imp
bio.imp = complete(imp.model,2)

saveRDS(bio.imp, file = paste0(save_folder,"bioImputedMICE.rds"))