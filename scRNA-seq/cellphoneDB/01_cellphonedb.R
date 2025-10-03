##/data/pingjiale/02_Software/01_anaconda/anaconda3/envs/R4.4.1/bin/Rscript
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
#library(SeuratData)
library(patchwork)
library(scater)
library(SeuratDisk)
set.seed(2)

blood <- readRDS("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/final_celltype.rds")
celltypes <- c("Tn_CD4", "Tcm_CD4","Th1_CD4","Th2_CD4","Th17_CD4","Treg","Tex_CD4",
               "Tn_CD8",'Tcm_CD8',"Tem","TCR","CTL","Tex_CD8",'proT',
               'NK1','NK2','Bn','Bm','Pla','cMC','nMC','cDC','pDC','HSPC')
Idents(blood) <- factor(blood$subtype, levels=celltypes)
donor <- unique(blood$donor)

for (i in donor){
        seu.sub <- subset(blood,  donor== i)
        rt=as.matrix(seu.sub@assays$RNA@data)
        gene=rownames(rt)
        rt2=cbind(gene,rt)

        dir.create(paste0('/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/07_celldb/',i))
        write.table(rt2, paste0('/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/07_celldb/',i,'/cellphonedb_count.txt'), sep='\t', quote=F,row.names = F,col.names = T)
        meta_data <- cbind(rownames(seu.sub@meta.data), seu.sub@meta.data[,'subtype', drop=F])
        meta_data <- as.matrix(meta_data)
        #meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
        write.table(meta_data, paste0('/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/07_celldb/',i,'/cellphonedb_meta.txt'), sep='\t', quote=F, row.names=F)
}

