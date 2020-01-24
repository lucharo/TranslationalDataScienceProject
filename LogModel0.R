
suppressPackageStartupMessages(library(tidyverse))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
covars = readRDS("data/Covars_toy.rds")

# Logistic model here

