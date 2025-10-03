#source activate R4.4.1
library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)

set.seed(2)
pbmc <- readRDS(file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/umap.rds')
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize")

cell.prop <- data.frame(table(pbmc$seurat_clusters,pbmc$donor))
colnames(cell.prop) <- c('cluster', 'donor', 'number')
cell.prop <- cell.prop %>% group_by(donor) %>% mutate(total=sum(number))
cell.prop$prop <- cell.prop$number/cell.prop$total
cell.prop <- cell.prop %>% group_by(cluster) %>% mutate(prop_total=sum(prop))
cell.prop$artio <-  cell.prop$prop / cell.prop$prop_total
cluster <- dcast(cell.prop, cluster~donor, value.var ='artio' )
write.csv(cluster,"/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/clusters_prop.csv", row.names = F)

p <- VlnPlot(pbmc, features = c("nFeature_RNA"),pt.size=0)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/nFeature.pdf", plot = p, width = 10, height = 4)
p <- VlnPlot(pbmc, features = c("percent.MT"),pt.size=0)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/mt.pdf", plot = p, width = 10, height = 4)

plan("multiprocess", workers = 6)
options(future.globals.maxSize = 100000 * 1024^5)
markers <- FindAllMarkers(pbmc, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = "wilcox")
write.csv(markers,file="/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/top_markers_all.csv",row.names=F)

VlnPlot <- VlnPlot(pbmc, features = c("PF4", "PPBP","TUBB1"),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/MEGA_V.png", plot = VlnPlot, width = 36, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c("PF4", "PPBP","TUBB1"), cols=c('grey90', '#C63C3C'),raster=FALSE,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/MEGA.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c("CD14", "FCGR3A",'MS4A7'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/MC_V.png", plot = VlnPlot, width = 36, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c("CD14", "FCGR3A",'MS4A7'), cols=c('grey90', '#C63C3C'),raster=FALSE,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/MC.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c("CMTM2", "CXCR2","FCGR3B","PTGS2",'CAMP','FCER1A'),pt.size = 0.011)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/Neu_V.png", plot = VlnPlot, width = 24, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c("CMTM2", "CXCR2","FCGR3B","PTGS2",'CAMP','FCER1A'), cols=c('grey90', '#C63C3C'),raster=FALSE,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/Neu.png", plot = FeaturePlot, width = 12, height = 18)

VlnPlot <- VlnPlot(pbmc, features = c("FCER1A", "CLEC10A","CD1C","IRF8"),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/DC_V.png", plot = VlnPlot, width = 24, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c("FCER1A", "CLEC10A","CD1C","IRF8"), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/DC.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c("FCGR3A", "NCAM1","KLRF1","NKG7"),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/NK_V.png", plot = VlnPlot, width = 24, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c("FCGR3A", "NCAM1","KLRF1","NKG7"), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/NK.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c("MS4A1", "CD19","MZB1","CD79A"),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/B_V.png", plot = VlnPlot, width = 24, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c("MS4A1", "CD19","MZB1","CD79A"), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/B.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c("CD27", "MZB1","CD38"),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/pla_V.png", plot = VlnPlot, width = 36, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c("CD27", "MZB1","CD38"), cols=c('grey90', '#C63C3C'),raster=FALSE,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/pla.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c('MME','IL4R','TCL1A','CLECL1','AIM2'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/BC_V.png", plot = VlnPlot, width = 36, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c('MME','IL4R','TCL1A','CLECL1','AIM2'), cols=c('grey90', '#C63C3C'),raster=FALSE,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/BC.png", plot = FeaturePlot, width = 12, height = 18)

VlnPlot <- VlnPlot(pbmc, features = c("CD3E", "CD3D","CD3G"),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/T_V.png", plot = VlnPlot, width = 36, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c("CD3E", "CD3D","CD3G"), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/T.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c('CD4','CD8A'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/TC_V.png", plot = VlnPlot, width = 24, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c('CD4','CD8A'), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/TC.png", plot = FeaturePlot, width = 12, height = 6)

VlnPlot <- VlnPlot(pbmc, features = c('CCR7','LEF1','AQP3'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/NaiveTTcmCD4+_V.png", plot = VlnPlot, width = 36, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c('CCR7','LEF1','AQP3'), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/NaiveTTcmCD4+.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c('CCR6','PRDM1','FOXP3','IL2RA'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/TemCD4+TregCD4+_V.png", plot = VlnPlot, width = 36, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c('CCR6','PRDM1','FOXP3','IL2RA'), cols=c('grey90', '#C63C3C'),raster=FALSE)#,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/TemCD4+TregCD4+.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c('GZMK','GZMB','GNLY','KLRB1','TRDV2'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/TemCD8+CTLCD8+_V.png", plot = VlnPlot, width = 36, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c('GZMK','GZMB','GNLY','KLRB1','TRDV2'), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/TemCD8+CTLCD8+.png", plot = FeaturePlot, width = 12, height = 18)

VlnPlot <- VlnPlot(pbmc, features = c('IFNG-AS1','GATA3','CCR6'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/Th_V.png", plot = VlnPlot, width = 36, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c('IFNG-AS1','GATA3','CCR6'), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/Th.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c('PDCD1','TIGIT','TOX'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/Exh_T_V.png", plot = VlnPlot, width = 36, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c('PDCD1','TIGIT','TOX'), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/Exh_T.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c('MKI67','TYMS'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/proTC_V.png", plot = VlnPlot, width = 24, height = 5)
FeaturePlot <- FeaturePlot(pbmc, features = c('MKI67','TYMS'), cols=c('grey90', '#C63C3C'),raster=FALSE,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/proTC.png", plot = FeaturePlot, width = 12, height = 6)

VlnPlot <- VlnPlot(pbmc, features = c('HBB','HBA1','HBA2','MALAT1'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/RBC_V.png", plot = VlnPlot, width = 36, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c('HBB','HBA1','HBA2','MALAT1'), cols=c('grey90', '#C63C3C'),raster=FALSE,order = TRUE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/RBC.png", plot = FeaturePlot, width = 12, height = 12)

VlnPlot <- VlnPlot(pbmc, features = c('KIT','PTPRC','CD34','PLXDC2','MEIS1','IGF2BP2'),pt.size = 0.01)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/HSPC_V.png", plot = VlnPlot, width = 36, height = 10)
FeaturePlot <- FeaturePlot(pbmc, features = c('KIT','PTPRC','CD34','PLXDC2','MEIS1','IGF2BP2'), cols=c('grey90', '#C63C3C'),raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/02_marker/gene_marker/HSPC.png", plot = FeaturePlot, width = 12, height = 18)
