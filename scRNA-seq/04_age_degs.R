#source activate R4.4.1
library(Seurat)
library(ggplot2)
library(MAST)
library(tidyverse)
#library(parallel)
set.seed(2)
library(Rcpp)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  NumericMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
' )

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#如果数据是浮点型需要把上述的IntegerMatrix替代为NumericMatrix，不然会强制转化为整型的矩阵。
pbmc <- readRDS(file = "/dellstorage02/quj_lab/pingjiale/02_Result/10_qz_pbmc/final_celltype.rds")
#pbmc <- readRDS("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/03_celltype/final_celltype.rds")
Idents(pbmc) <- pbmc$subtype

pbmc <- subset(pbmc, group  == 'control')
celltypes <- c("Tn_CD4", "Tn_CD8","Th2_CD4",'Tcm_CD8',"Tem","TCR","CTL","Tex_CD8",'proT',
			   "Tcm_CD4","Th1_CD4","Th17_CD4","Treg","Tex_CD4",
			   'NK1','NK2','Bn','Bm','Pla','cMC','nMC','cDC','pDC','HSPC')

df <- data.frame()

# 使用mclapply并行化循环
for (i in celltypes) {
  tmp <- subset(pbmc, subtype == i)
  #Idents(tmp) <- tmp$Age
  #set.seed(2024)
  #tmp <- subset(x = tmp, downsample = 1000) #no downsample —— Killed
  
  Idents(tmp) <- tmp$donor
  set.seed(2024)
  tmp <- subset(x = tmp, downsample = 400) #no downsample —— Killed
  Idents(tmp) <- tmp$subtype

  age_num = as.numeric(tmp$Age)
  log_counts <- tmp@assays$RNA@data

  #log_counts[log_counts < 0.1] <- 0
  log_counts <- as(log_counts, "dgCMatrix")
  s = as_matrix(log_counts)
  s = s[which(rowSums(s) > 0),]

  fData = data.frame(names=rownames(s))
  rownames(fData) = rownames(s)

  cData = data.frame(cond=age_num,gender=tmp$gender)
  rownames(cData) = colnames(log_counts)
  
  obj <- FromMatrix(s, cData, fData)#,check_sanity = FALSE)

  colData(obj)$cngeneson <- scale(colSums(assay(obj)>0)) #对每个细胞中检测到的基因总数进行矫正
  zlmCond <- zlm(formula = as.formula("~cond+gender+cngeneson"), obj, parallel = TRUE)

  summaryCond <- summary(zlmCond, doLRT="cond")
  summaryDt <- summaryCond$datatable
  dt1 = summaryDt[contrast=="cond" & component=="H", .(primerid, `Pr(>Chisq)`)]
  dt2 = summaryDt[contrast=="cond" & component=="logFC", .(primerid, coef, z)]
  de_res = merge(dt1, dt2, by="primerid")
  colnames(de_res) <- c("gene", "age.H_p", "age.logFC_coef", 'age.logFC_z')
  de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")

  de_res$subtype <- i
  write.csv(de_res, paste0("/dellstorage02/quj_lab/pingjiale/02_Result/10_qz_pbmc/age_degs_gender/",i,"_gene_downsample_donor.csv"))
  df <- rbind(df, de_res)
  print(paste0(i, ' is finished'))
}
# 合并结果
df_deg <- subset(df, age.H_fdr < 0.05 & abs(age.logFC) > 0)
write.csv(df, "/dellstorage02/quj_lab/pingjiale/02_Result/10_qz_pbmc/age_degs_gender/all_age_gene.csv")
write.csv(df_deg, "/dellstorage02/quj_lab/pingjiale/02_Result/10_qz_pbmc/age_degs_gender/related_age_degs.csv")
