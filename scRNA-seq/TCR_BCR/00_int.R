#source activate R4.4.1
library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

set.seed(2)
bc <- readRDS("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/08_TBcr/female.combineB.rds")
bc$barcode <- colnames(bc)
bcr <- bc@meta.data[,c("barcode",'batch',"donor","gender","group","Age","subtype","CTgene","CTaa")]
bcr_female <- subset(bcr, ! CTaa == "NA")

bc <- readRDS("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/08_TBcr/male.combineB.rds")
bc$barcode <- colnames(bc)
bcr <- bc@meta.data[,c("barcode",'batch',"donor","gender","group","Age","subtype","CTgene","CTaa")]
bcr_male <- subset(bcr, ! CTaa == "NA")
bcr <- rbind(bcr_female,bcr_male)
write.csv(bcr,"/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/08_TBcr/bcr_meta.csv", row.names = F)

tc <- readRDS("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/08_TBcr/female.combineT.rds")
tc$barcode <- colnames(tc)
tcr <- tc@meta.data[,c("barcode",'batch',"donor","gender","group","Age","subtype","CTgene","CTaa")]
tcr_female <- subset(tcr, ! CTaa == "NA")

tc <- readRDS("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/08_TBcr/male.combineT.rds")
tc$barcode <- colnames(tc)
tcr <- tc@meta.data[,c("barcode",'batch',"donor","gender","group","Age","subtype","CTgene","CTaa")]
tcr_male <- subset(tcr, ! CTaa == "NA")
tcr <- rbind(tcr_female,tcr_male)
write.csv(tcr,"/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/08_TBcr/tcr_meta.csv", row.names = F)
