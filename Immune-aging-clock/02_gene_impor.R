#install.packages('glmnet')
library(glmnet)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(cowplot)
library(hstats)
library(foreach)
library(doParallel)

########train or valid#####
inf <- read.csv("E:/01_progrom/06_qz_pbmc_new/00_inf/inf_merge_final_me.csv")
celltype <- read.csv("E:/01_progrom/06_qz_pbmc_new/03_celltype/final_celltype_prop.csv",row.names = 1)
inf <- subset(inf, donor %in% unique(celltype$donor))
inf <- subset(inf, group  == 'control')

inf$stage <- cut(inf$Age, 
                 breaks = c(-Inf, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, Inf), 
                 labels = c("stage_01", "stage_02", "stage_03", 
                            "stage_04", "stage_05", "stage_06",
                            "stage_07", "stage_08", "stage_09", 
                            "stage_10", "stage_11", "stage_12", "stage_13"))

set.seed(2024)
sampled_prop <- inf %>% group_by(gender,stage) %>% sample_frac(size = 0.5)
#aa <- c('male-pbmc-28','male-pbmc-15','male-pbmc-102','206868560011_R02C01')
#bb <- intersect(sampled_prop$donor,aa)

valid.rows <- sampled_prop$donor
train.rows <- setdiff(inf$donor, sampled_prop$donor)

########module#####
inf <- read.csv("E:/01_progrom/06_qz_pbmc_new/00_inf/inf_merge_final_me.csv")
inf <- inf[,c('donor', 'Age','gender')]

celltypes <- c("Tn_CD4", "Tcm_CD4","Th1_CD4","Th2_CD4","Th17_CD4","Treg",
               "Tex_CD4","Tn_CD8",'Tcm_CD8',"Tem","TCR","CTL",
               "Tex_CD8",'proT','NK1','NK2','Bn','Bm',
               'Pla','cMC','nMC','cDC','pDC','HSPC',
               'CD4T','CD8T')
annogene <- read.csv('E:/01_progrom/06_qz_pbmc_new/001_analysis/03_clock_ref/hg38_Gene_info_with_Exon_len.csv')
annogene <- subset(annogene, ! (gene_biotype %in% c('artifact','misc_RNA')))

registerDoParallel(cores = detectCores() - 1)
packs <- c("Hmisc",'ggplot2','glmnet','dplyr','ggthemes','limma','cowplot','hstats')
foreach(cell = celltypes, .packages = packs) %dopar% {
  #for (cell in celltypes){
  meta <- read.csv(paste0("E:/01_progrom/06_qz_pbmc_new/05_gene/expression_donor/cor_age_gender_sample_expression/",cell,"_cor_gene_meta.csv"),check.names = FALSE)
  aa <- grep("^HIST", setdiff(colnames(meta), annogene$gene_name), value = TRUE, invert = TRUE)
  meta <- meta[, -which(colnames(meta) %in% aa[5:length(aa)])]
  
  valid.rows <- intersect(valid.rows,meta$donor)
  train.rows <- intersect(train.rows,meta$donor)
  
  wtmeta <- subset(meta, group %in% "control")
  hrtmeta <- subset(meta, group  == 'Used')
  hrtmeta = hrtmeta[,c(1,4:length(colnames(hrtmeta)))]
  rownames(hrtmeta) <- hrtmeta$donor
  hrtmeta <- hrtmeta[,-1]
  
  mae.df=c()
  dir.create(paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/02_gene/',cell,'/'))
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/02_gene/',cell,'/')

  mae <- read.csv(paste0(result,'mae_result.csv'))
  mae <- mae[order(mae$MAE),]
  alp <- mae$alpha[1]

    set.seed(2024)
    #dir.create(paste0(result,alp))
    pd <- unique(wtmeta[,c('donor','Age', 'gender')])
    pd$dataset <- pd$donor
    pd[pd$dataset %in% train.rows,]$dataset <- 'Train'
    pd[pd$dataset %in% valid.rows,]$dataset <- 'Valid'
    
    shannodata= wtmeta[,c(1,4:length(colnames(wtmeta)))]
    rownames(shannodata) <- shannodata$donor
    shannodata <- shannodata[,-1]##female=0,male=1
    x=shannodata[train.rows,]
    
    y = subset(pd, donor %in% train.rows)
    rownames(y) <- y$donor
    y <- y[rownames(x),]
    y = as.numeric(y$Age)
    
    fitcv <- cv.glmnet(as.matrix(x),y,alpha= alp)
    model = glmnet(as.matrix(x),y, family="gaussian", alpha=alp, lambda = fitcv$lambda.min)
    
    s <- perm_importance(model, y = y, X = as.matrix(x))
    plot(s[1:10])+theme_bw()
    ggsave(paste0(result,alp,'/perm_importance.pdf'),width =2.5,height = 2.5)
    write.csv(s$M,paste0(result,alp,'/importance.csv'))
    write.csv(s$SE,paste0(result,alp,'/importance_SE.csv'))
}
