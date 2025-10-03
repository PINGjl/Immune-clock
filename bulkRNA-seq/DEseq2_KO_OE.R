library(DESeq2)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)

####all_pac#####
data.path="E:/01_progrom/06_qz_pbmc_new/10_bulk/01_count/CD4_OE/"
sample.path="E:/01_progrom/06_qz_pbmc_new/10_bulk/02_result/CD4_OE/"
data.list=list.files(data.path,pattern = '*txt',full.names = T)
data.list
##确定字符串数
##df=read.table(data.list[1],header = F)
##name=substr(data.list[1],43,nchar(data.list[1])-8)
##head(name)   "KOP15_1"
##确定去掉最后几行(乱码)
##df=df[-((nrow(df)-4):nrow(df)),]  去掉的是最后5行
##tail(df,10)
HTseq.handling=function(x,n){
  df=read.table(x,header = F)
  name=substr(x,n,nchar(x)-8)
  colnames(df)=c("ensembl_gene_id",name)
  df=df[-((nrow(df)-4):nrow(df)),]
  return(df)
}
tt=lapply(data.list, HTseq.handling,54)
head(tt[[1]])
tail(tt[[1]])
head(tt[[6]])
data.name=unique(substr(data.list,54,nchar(data.list)-10))
allmerge=function(x){
  merge.data=merge(x[[1]],x[[2]],by='ensembl_gene_id')
  for (i in 3:length(x)) {
    merge.data=merge(merge.data,x[[i]],by='ensembl_gene_id')
  }
  return(merge.data)
}
data.all=allmerge(tt)

data.combn=data.frame(c(data.name[1],data.name[2]))
all.name=as.character(data.combn[,1])
all.sample=c(paste(all.name[1],c(1:3),sep='-'),paste(all.name[2],c(1:3),sep='-'))

data.all=data.all[,c("ensembl_gene_id", all.sample)]##换列
write.csv(data.all,paste(sample.path,'merge_sample.csv',sep = ''),row.names =F)

samplecondition=factor(c(rep(all.name[1],3),rep(all.name[2],3)),
                       levels = c(all.name[1],all.name[2]))
samplecondition ##差异比较矩阵,说明是对照还是处理
mycol.data=data.frame(row.names = all.sample,condition=samplecondition)
mycol.data ##样品信息矩阵，即condition
mycount.data=data.all
rownames(mycount.data)=mycount.data[,1]
mycount.data=mycount.data[,-1] ##表达矩阵（所有样品）
class(mycount.data)

##构建dds矩阵
mydds=DESeqDataSetFromMatrix(countData = mycount.data,colData = mycol.data,design = ~ condition)
head(mycount.data)
head(mycol.data)
samplecondition

##标准化
mydds=DESeq(mydds)
myres=results(mydds)
########rlog
rld=rlog(mydds)
write.csv(assay(rld),paste(sample.path,"merge_sample_rlog.csv",sep=''),row.names = T)
myres=myres[order(myres$padj),]

##画图
print('Euclidean distance plot')
head(rld)
head(assay(rld))
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- all.sample
rownames(sampleDistMatrix)
colnames(sampleDistMatrix) <-rownames(sampleDistMatrix)
sampleDist_colors <- colorRampPalette(rev(brewer.pal(9, 'Blues')))(200)
pdf(paste(sample.path,"Euclidean_distances_heatmap.pdf",sep=""),width=6, height=4)
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           cluster_rows = T,
           cluster_col = T,
           treeheight_col = 0,
           treeheight_row = 100,
           col=sampleDist_colors,
           #breaks = seq(0,300,length.out = length(sampleDist_colors)),
           border_color=NA,
           cellwidth = 10,
           cellheight = 10)
print(p)
dev.off()

print('principal components analysis and plot')
pca_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
write.csv(pca_data,paste(sample.path,"PCA_data.csv",sep=""),row.names = F)
#pca_data <- read.csv(paste(sample.path,"PCA_data.csv",sep=""))
percentVar <- round(100 * attr(pca_data, "percentVar"))
pdf(paste(sample.path,"PCA_plot.pdf",sep=""),width=4, height=2.6)
p=ggplot(pca_data, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) + 
  #scale_color_manual(values = c('#CEB7B3', '#47B1B6', '#E6949A'))+ 
  scale_shape_manual(values = c(16,16,16,16))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+theme_classic()#+
  #stat_ellipse(level = 0.95)#+stat_ellipse(level = 0.8)
print(p)
dev.off()


all_symbol=read.csv("E:/01_progrom/04_lung_TB/02_bulkRNA/Homo_sapiens.GRCh37.87_chr_with_len.tsv",sep='\t')
##提取差异分析结果
all_gene.data=merge(as.data.frame(myres),as.data.frame(counts(mydds,normalize=TRUE)),by="row.names",sort=FALSE)
resdata=all_gene.data
colnames(all_symbol)[1]="gene_id"
colnames(resdata)[1]="gene_id"
resdata_anno=merge(all_symbol,resdata,by="gene_id")
write.csv(resdata_anno,paste(sample.path,"all_gene.csv",sep=''),row.names = F)
######diff_gene
diff_gene=subset(resdata_anno,padj<0.05 & (log2FoldChange>=1|log2FoldChange<=-1)) #padj
write.csv(diff_gene,paste(sample.path,"diff_gene.csv",sep=''),row.names = F)
#####up_gene
up_gene=subset(diff_gene,log2FoldChange>=1)
write.csv(up_gene,paste(sample.path,"up_gene.csv",sep=''),row.names = F)
#####down_gene
down_gene=subset(diff_gene,log2FoldChange<=-1)
write.csv(down_gene,paste(sample.path,'down_gene.csv',sep=''),row.names = F)


####cut tag###
all_symbol=read.csv("E:/01_progrom/04_lung_TB/02_bulkRNA/Homo_sapiens.GRCh37.87_chr_with_len.tsv",sep='\t')
##提取差异分析结果
peak <- read.csv('E:/01_progrom/06_qz_pbmc_new/11_CutTag/Up_peak_anno.csv')

common <- merge(all_symbol,peak,by.x = "GeneID",by.y = "geneId")
write.csv(common,'E:/01_progrom/06_qz_pbmc_new/11_CutTag/Up_peak_anno_gene.csv',row.names = F)
