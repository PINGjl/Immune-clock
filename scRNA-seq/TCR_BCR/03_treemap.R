library(ggplot2)
library(treemap)
library(randomcoloR)
library(RColorBrewer)
library(groupdata2)
#####TC and BC#####
data <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/bcr_meta.csv')
data <- subset(data, group == 'control')
#data <- subset(data, donor == 'male-pbmc-24')
data$group <- as.factor(ifelse(data$Age <= 29,'stage_1',
                                  ifelse(data$Age <= 39 ,'stage_2',
                                      ifelse(data$Age <= 49 ,'stage_3',
                                          ifelse(data$Age <= 59 ,'stage_4',
                                              ifelse(data$Age <= 69 ,'stage_5','stage_6'))))))

data$group <- as.factor(ifelse(data$Age <= 39,'stage_1',
                               ifelse(data$Age <= 59 ,'stage_2','stage_3')))

sample <- unique(data[,c('group','donor')])
sample <- subset(sample, ! donor %in% c('male-pbmc-24','male-pbmc-27','206863650120_R02C02'))
set.seed(2025)
samples <- sample %>% group_by(group) %>% slice_sample(n = 20)
data <- subset(data, donor %in% samples$donor)

###感觉不能只看克隆扩增的
#freq <- as.data.frame(table(data$donor,data$CTaa))
#freq <- subset(freq, Freq > 1)
#colnames(freq) <- c('donor','CTaa','count')
#data <- merge(data, freq, by.x = c('donor','CTaa'), by.y = c('donor','CTaa'))

#set.seed(2025)
#data <- downsample(data, cat_col = 'group')
num <- as.data.frame(table(data$donor,data$CTaa))

data <- merge(data, num, by.x = c('donor','CTaa'), by.y = c('Var1','Var2'))
#stage_1 stage_2 stage_3 
#22065   20238   21467
#stage_1 stage_2 stage_3 
#5372    6091    5247
for (i in c('stage_1','stage_2','stage_3')) { #,'stage_4','stage_5','stage_6'
  tmp <- subset(data, group == i)
  tmp <- unique(tmp[,c('donor','CTaa','Freq')])
  tmp <- tmp[order(tmp$Freq,decreasing = T),]
  
  tmp$ident <- sample(c(1:max(length(tmp$donor))))
  palette <- distinctColorPalette(60) #差异明显的60种
  pdf(paste0("E:/01_progrom/06_qz_pbmc_new/08_TBcr/03_treemap/bc_",i,".pdf"),width=5, height=5)
  set.seed(2025)
  treemap(tmp, index=c('donor','CTaa'), vSize="Freq", 
        vColor="ident",
        type="manual",palette = palette, #manual
        border.col = c("black","white"),#边框颜色
        border.lwds = 0.2,#边框线宽度
        title=i)
  dev.off()
}
