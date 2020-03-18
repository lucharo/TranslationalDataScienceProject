#Script to correct patient IDs in SNP datasets (both snpProcessed and snpImputed), 
#using Barbz bridge file

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cluster = 1

if (cluster == 1){
  save_data = data_folder = "../FULLDATA/preprocessed/"
  save_plots = "../FULLResults/"
} else {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  save_data = data_folder = "../data/preprocessed/"
  save_plots = "../results/"
}

snp = readRDS(paste0(data_folder,"snpProcessed.rds"))
snpImp = readRDS(paste0(data_folder,"snpImputed.rds"))
bf = readRDS(paste0(data_folder,"bridge_file.rds"))

bf = as.data.frame(bf)

snp.new = merge(bf, snp, by.x='previous', by.y='ID')
snp.new = snp.new[,-1]
colnames(snp.new)[1] = 'ID'

snpImp.new = merge(bf, snpImp, by.x='previous', by.y='ID')
snpImp.new = snpImp.new[,-1]
colnames(snpImp.new)[1] = 'ID'

saveRDS(snp.new, paste0(save_data,"snpProcessed.rds"))
saveRDS(snpImp.new, paste0(save_data,"snpImputed.rds"))


