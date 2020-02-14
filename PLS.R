rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

bio <- readRDS("data/preprocessed/bioImputed.rds")
