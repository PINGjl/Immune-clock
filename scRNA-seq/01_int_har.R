#source activate R4.4.1
library(Seurat)
library(patchwork)
library(cowplot)
library(harmony)
library(ggplot2)
library(reshape2)
library(dplyr)
library(future)

data <- readRDS(file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/1_merge_male_female.rds") #merge.raw.RDS
data[["percent.MT"]] <- PercentageFeatureSet(data, pattern = "^MT-")
num <- as.data.frame(table(data$donor,data$gender))
num <- subset(num, Freq > 600)
data <- subset(data, donor %in% num$Var1)
p <- VlnPlot(data, features = c("nFeature_RNA"), group.by = 'donor', ncol = 1)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/gene_vln.pdf", plot = p,limitsize = FALSE, width = 50, height = 10)
p <- VlnPlot(data, features = c("percent.MT"), group.by = 'donor', ncol = 1)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/mt_vln.pdf", plot = p,limitsize = FALSE, width = 50, height = 10)

data_qc <- subset(data,subset = percent.MT < 5 & nFeature_RNA > 400 & nFeature_RNA < 4000) 

inf <- read.csv('/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/inf_merge_final_me.csv')
obj<-inf$Age
names(obj) <- inf$donor
data_qc@meta.data <- data_qc@meta.data %>% mutate(Age=obj[donor])
obj<-inf$group
names(obj) <- inf$donor
data_qc@meta.data <- data_qc@meta.data %>% mutate(group=obj[donor])

#####质控结果#####
zhikong <- data_qc@meta.data[,c("donor",'Age','gender','group',"nFeature_RNA", "nCount_RNA")]
zhikong <- zhikong %>% group_by(donor) %>% mutate(gene_num=mean(nFeature_RNA))
zhikong <- zhikong %>% group_by(donor) %>% mutate(UMI_num=mean(nCount_RNA))
zhikong <- unique(zhikong[,c("donor",'Age','gender','group',"gene_num", "UMI_num")])
number <- as.data.frame(table(data_qc$donor))
colnames(number) <- c("donor",'cell_num')
zhikong <- merge(zhikong, number, by = "donor")
write.csv(zhikong,"/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/qc_per_donor.csv", row.names = FALSE)

#####one way
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^5)
data_qc <- SCTransform(data_qc, verbose = FALSE) ###这种方法超级慢 SCTransform=Normalization+FindVariableFeatures+ScaleData 会好?
data_qc <- RunPCA(data_qc, verbose=FALSE,npcs = 50,features =VariableFeatures(data_qc))
saveRDS(data_qc, file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/SCT/SCT_Data.rds')
data_har <- data_qc %>% RunHarmony("batch", assay.use = "SCT") #Harmony converged after 6 iterations

#####two way -final
data_qc <- NormalizeData(data_qc)
data_qc <- FindVariableFeatures(data_qc, selection.method = "vst",nfeatures = 3000)#默认2000
data_qc <- ScaleData(data_qc)
data_qc <- data_qc %>% RunPCA(npcs = 50, verbose = FALSE,features =VariableFeatures(data_qc))
saveRDS(data_qc, file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/Normalize_Data.rds')
data_har <- data_qc %>% RunHarmony(group.by.vars = "batch") #Harmony converged after 5 iterations

p <- DimPlot(object = data_qc, reduction = "pca", pt.size = .1, group.by = "batch",raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/harmony_batch_be.pdf", plot = p, width = 6, height = 5)
p <- DimPlot(object = data_har, reduction = "harmony", pt.size = .1, group.by = "batch",raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/harmony_batch.pdf", plot = p, width = 6, height = 5)

pdf('/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/pca_heatmap.pdf', height=20, width=10)
DimHeatmap(data_har, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()
pdf('/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/pca_elbow.pdf', height=5, width=8)
p=ElbowPlot(data_har, ndims = 50)
print(p)
dev.off()
saveRDS(data_har, file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/before_pc.rds')

pc.num=c(1:20)
#pc.num=c(1:20)
data_har <- readRDS(file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/before_pc.rds')

data_har <- FindNeighbors(data_har,reduction = "harmony", dims = pc.num)
data_har <- FindClusters(data_har, resolution = 2.0)#2.0
data_har <- RunUMAP(data_har,reduction = "harmony", dims = pc.num)#,spread = 3)

plot <- DimPlot(data_har, reduction = "umap", label = TRUE,pt.size = 0.1,raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/umap.pdf", plot = plot, width = 8, height = 6)
plot <- DimPlot(data_har, reduction = "umap", group.by = "donor",label = TRUE,pt.size = 0.1,raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/umap_donor.pdf", plot = plot, width = 28, height = 5)
saveRDS(data_har, file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/umap.rds')

#SCT —— 效果不好
pc.num=c(1:15,20)
data_har <- readRDS(file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/SCT/before_pc.rds')

data_har <- FindNeighbors(data_har,reduction = "harmony", dims = pc.num)
data_har <- FindClusters(data_har, resolution = 2.0)#2.0
data_har <- RunUMAP(data_har,reduction = "harmony", dims = pc.num)#,spread = 3)

plot <- DimPlot(data_har, reduction = "umap", label = TRUE,pt.size = 0.1,raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/SCT/umap.pdf", plot = plot, width = 8, height = 6)

plot <- DimPlot(data_har, reduction = "umap", group.by = "donor",label = TRUE,pt.size = 0.1,raster=FALSE)
ggsave("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/SCT/umap_donor.pdf", plot = plot, width = 28, height = 5)

plot <- DimPlot(data_har, reduction = "umap", split.by = "donor",pt.size = 0.1,label.size = 6,ncol=16,raster=FALSE)
ggsave('/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/SCT/umap_sample_single.pdf',plot = plot,width=30,height=30)
saveRDS(data_har, file = '/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/01_int/SCT/umap.rds')
