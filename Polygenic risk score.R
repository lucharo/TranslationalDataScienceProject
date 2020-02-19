# Aim of this script is to calculate a polygenic risk score for each person

##################################################################
##                 Prepare libraries and data                   ##
##################################################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

snp_info <- readRDS("data/preprocessed/snpInfo.rds")
snp <- readRDS("data/preprocessed/snpProcessed.rds")
