#source activate R4.4.1
library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)

set.seed(2)
pbmc <- readRDS(file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/umap.rds')
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize")

pbmc <- subset(pbmc, idents = c(41,42,53,54,43,49,48,52), invert = TRUE)

#Idents(pbmc) <- pbmc$seurat_clusters
celltypes <- c('CD4T','CD8T','NK1','NK2','Bn','Bm','Pla','cMC','nMC','cDC','pDC','HSPC')
C1 <- c(0,2,4,6,8,12,16,19,21,22,24,25,32,36,37,38)
C2 <- c(1,3,5,10,13,15,17,18,20,26,44,46,47,51)
C3 <- c(31)
C4 <- c(7,9,14,28)

C5 <- c(11)
C6 <- c(23,27,30)
C7 <- c(40)

C8 <- c(29)
C9 <- c(33,34)
C10<- c(39)
C11<- c(45)
C12 <- c(35,50)

for (i in 1:12) {
  m <- paste0('C',i)
  tmp <- get(m)
  assign(celltypes[i],tmp)
}

new.cluster.ids <- c(rep(NA,length(unique(Idents(pbmc))  )))
cluster <- function(celltypes,new.cluster.ids){
  for (cell in celltypes) {
    tmp <- get(cell)
    for (i in tmp) {
      new.cluster.ids[i+1] <- cell
    }
  }
  return(new.cluster.ids)

}
new.cluster.ids <- cluster(celltypes,new.cluster.ids)

pbmc <- RenameIdents(pbmc,c("0"=new.cluster.ids[1], "1"=new.cluster.ids[2], "2"=new.cluster.ids[3], "3"=new.cluster.ids[4],"4"= new.cluster.ids[5], "5"=new.cluster.ids[6],"6"=new.cluster.ids[7], "7"=new.cluster.ids[8],"8"= new.cluster.ids[9],"9"= new.cluster.ids[10],"10"= new.cluster.ids[11],"11"=new.cluster.ids[12], "12"=new.cluster.ids[13], "13"=new.cluster.ids[14],"14"= new.cluster.ids[15],"15"=new.cluster.ids[16],"16"=new.cluster.ids[17], "17"=new.cluster.ids[18],"18"= new.cluster.ids[19],"19"= new.cluster.ids[20],"20"= new.cluster.ids[21],"21"=new.cluster.ids[22], "22"=new.cluster.ids[23], "23"=new.cluster.ids[24],"24"= new.cluster.ids[25],"25"=new.cluster.ids[26],"26"=new.cluster.ids[27], "27"=new.cluster.ids[28],"28"= new.cluster.ids[29],"29"= new.cluster.ids[30],"30"= new.cluster.ids[31],"31"=new.cluster.ids[32], "32"=new.cluster.ids[33], "33"=new.cluster.ids[34],"34"= new.cluster.ids[35],"35"=new.cluster.ids[36],"36"=new.cluster.ids[37], "37"=new.cluster.ids[38],"38"= new.cluster.ids[39],"39"= new.cluster.ids[40],"40"= new.cluster.ids[41],"41"=new.cluster.ids[42], "42"=new.cluster.ids[43], "43"=new.cluster.ids[44],"44"= new.cluster.ids[45], "45"=new.cluster.ids[46],"46"=new.cluster.ids[47], "47"=new.cluster.ids[48],"48"= new.cluster.ids[49],"49"= new.cluster.ids[50],"50"= new.cluster.ids[51],"51"= new.cluster.ids[52]))

Idents(pbmc) <- factor(Idents(pbmc), levels=celltypes)
pbmc$celltype <- Idents(pbmc)

CT.col <- c("#75BC57","#2DB7B8","#45d9fd","#84B1ED",'#D499B9','#9055A2','#5c196b','#F17F42','#ED8272','#f9c00c','#ddd954',"#f1404b")#"#4ea1d3",

UMAP <- DimPlot(pbmc, reduction = "umap", cols=CT.col, label = TRUE, pt.size = 0.01,label.size = 3.0,raster = FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/celltype_UMAP.pdf", plot = UMAP, width = 6, height = 5)
UMAP <- DimPlot(pbmc, reduction = "umap", cols=CT.col, label = FALSE, pt.size = 0.01,raster = FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/celltype_UMAP-2.pdf", plot = UMAP, width = 6, height = 5)

#'IGF2BP2','MKI67','TYMS',
markers <- c('CD3D','CD3E','CD3G','CD4','CD8A','KLRF1','NCAM1','FCGR3A','MS4A1','CD79A','TCL1A','IL4R','CLECL1','AIM2','MZB1','CD38','MS4A7','CD14','FCER1A','CLEC10A','CD1C','IRF8','MEIS1','CD34')
Idents(pbmc) <- factor(pbmc$celltype, levels=rev(celltypes))
DotPlot<-DotPlot(pbmc,features= markers, cols=c('grey90', '#C63C3C'))+RotatedAxis()
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/DotPlot.pdf", plot = DotPlot, width = 10, height = 4)

saveRDS(pbmc, file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/Combination_celltype.rds")

###prop####
cell_number=table(pbmc$celltype,pbmc$donor)
cell.prop <- data.frame(table(pbmc$celltype,pbmc$donor))
colnames(cell.prop) <- c('celltype', 'donor', 'number')
cell.prop <- cell.prop %>% group_by(donor) %>% mutate(total=sum(number))
cell.prop$prop <- cell.prop$number/cell.prop$total
write.csv(cell.prop,"/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/celltype_prop.csv")

####anno合并###
all.meta <- readRDS(file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/Combination_celltype.rds")
all.meta@meta.data$subtype <- rownames(all.meta@meta.data)

cd4.meta <- readRDS(file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/cd4tc_har/Combination_celltype.rds")
info <- as.data.frame(rownames(cd4.meta@meta.data))
info$subtype <- cd4.meta@meta.data$celltype
numbers <- length(info$subtype)
for (i in c(1:numbers)) {
  p <- info[i,1]
  g <- info[i,2]
  all.meta@meta.data$subtype[which(all.meta@meta.data$subtype == p, arr.ind=T)] <- paste(g)
}

cd8.meta <- readRDS(file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/cd8tc_har/Combination_celltype.rds")
info <- as.data.frame(rownames(cd8.meta@meta.data))
info$subtype <- cd8.meta@meta.data$celltype
numbers <- length(info$subtype)
for (i in c(1:numbers)) {
  p <- info[i,1]
  g <- info[i,2]
  all.meta@meta.data$subtype[which(all.meta@meta.data$subtype == p, arr.ind=T)] <- paste(g)
}

sub.meta <- subset(all.meta, celltype %in% c('NK1','NK2','Bn','Bm','Pla','cMC','nMC','cDC','pDC','HSPC'))
info <- as.data.frame(rownames(sub.meta@meta.data))
info$subtype <- sub.meta@meta.data$celltype
numbers <- length(info$subtype)
for (i in c(1:numbers)) {
  p <- info[i,1]
  g <- info[i,2]
  all.meta@meta.data$subtype[which(all.meta@meta.data$subtype == p, arr.ind=T)] <- paste(g)
}

fin.meta <- subset(all.meta, subtype %in% c("Tn_CD4", "Tcm_CD4","Th1_CD4","Th2_CD4","Th17_CD4","Treg","Tex_CD4",
										"Tn_CD8",'Tcm_CD8',"Tem","TCR","CTL","Tex_CD8",'proT',
										'NK1','NK2','Bn','Bm','Pla','cMC','nMC','cDC','pDC','HSPC'))
Idents(fin.meta) <- fin.meta$subtype
saveRDS(fin.meta, file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/final_celltype.rds")

UMAP <- DimPlot(fin.meta, reduction = "umap", label = TRUE, pt.size = 0.1,label.size = 3.0,raster = FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/final_celltype.pdf", plot = UMAP, width = 8, height = 5)
###prop####
cell.prop <- data.frame(table(fin.meta$subtype,fin.meta$donor))
colnames(cell.prop) <- c('subtype', 'donor', 'number')
cell.prop <- cell.prop %>% group_by(donor) %>% mutate(total=sum(number))
cell.prop$prop <- cell.prop$number/cell.prop$total
write.csv(cell.prop,"/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/final_celltype_prop.csv")
