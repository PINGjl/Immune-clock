###pace####
get_age_pace = function(CA, predAge){
  lm = lm(predAge~CA)
  pace = predAge - (CA * lm$coefficients[2] + lm$coefficients[1])
  return(pace)
}

pred.all <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/predAge_all.csv')

df <- c()
type <- unique(pred.all$subtype)
for (i in type){
  pred <- subset(pred.all, subtype == i)
  
  pace <- data.frame(sample = pred$donor, Age = pred$Age, gender = pred$gender, pred = pred$pred_age, 
                     Age_pace = get_age_pace(pred$Age,pred$pred_age))
  pace$subtype = i
  df <- rbind(df, pace)
  
  p = ggplot() +
    geom_point(data = pace, aes(Age, Age_pace), size = 1) +
    #geom_abline(data = NULL,linetype=4,color="grey") +
    geom_smooth(data = pace, aes(Age, Age_pace), method = lm, se = F, color = '#3d60b7') +
    #scale_y_continuous(limits = c(10, 75)) +
    #scale_x_continuous(limits = c(10, 95)) +
    labs(y = paste0('Age_pace')) +
    labs(x = paste0('Age')) +
    theme_base() +
    theme(legend.position = 'none')
  ggsave(paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/plot/',i,'_val_age_result_pace.pdf'), plot = p, width = 4, height =4)#7
}
write.csv(df,'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/Age_pace_all.csv', row.names = F)

#####age 加速和减速 top20% 聚类#####
pace_value <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/Age_pace_all.csv')
pace_value <- subset(pace_value, !subtype %in% c('CD4T','CD8T',"all")) #"all""all_sample""TCRAge"
type <- unique(pace_value$subtype)

df <- c()
for (i in type){
  pace <- subset(pace_value, subtype == i)
  pace <- pace[order(pace$Age_pace),]
  up <- pace$sample[c((length(pace$sample)-20):length(pace$sample))]
  down <- pace$sample[c(1:21)]
  pace$group[pace$sample %in% up] <- 1
  pace$group[pace$sample %in% down] <- -1
  pace <- subset(pace, !group == '')
  df <- rbind(df, pace)
}

####num
num <- as.data.frame(table(df$sample,df$group))
num <- dcast(num, Var1 ~ Var2, value.var = "Freq")
num$Freq <- num$`1` - num$`-1`
num <- subset(num, Freq >= 7 | Freq <= -7)
num$Group <- as.factor(ifelse(num$Freq > 0 ,'Age-accelerators','Age-decelerators'))
table(num$Group)

###Data for Heatmap
aa.pace=dcast(df, sample ~ subtype, value.var = "Age_pace")
rownames(aa.pace)=aa.pace$sample
aa.pace=aa.pace[,-1]
nrow(aa.pace)-colSums(is.na(aa.pace))
aa <- setdiff(type,c("pAge","TCRAge","tc","nk","bc","mc","dc","all_sample"))
order.list=c("pAge",aa,"TCRAge","tc","nk","bc","mc","dc","all_sample")

aa.pace=aa.pace[,order.list]
aa.pace$diff=rowSums(aa.pace>0,na.rm = T) - rowSums(aa.pace<0,na.rm = T)
###
n=32
up.pace=aa.pace[aa.pace$diff>=0,]
up.pace=up.pace[order(up.pace$diff,rowSums(up.pace[,-c(n+1)]),decreasing = T),]
up.pace$order1=apply(up.pace[,-n],1, function(x) min(which(x>0)))
up.pace$order2=apply(up.pace[,-n],1, function(x) sum(which(x>0)))
up.pace=up.pace[order(up.pace$diff,-up.pace$order1,-up.pace$order2,decreasing = T),]
###
down.pace=aa.pace[aa.pace$diff<0,]
down.pace=down.pace[order(down.pace$diff,rowSums(down.pace[,-c(n+1)]),decreasing = T),]
down.pace$order1=apply(down.pace[,-n],1, function(x) min(which(x<0)))
down.pace$order2=apply(down.pace[,-n],1, function(x) sum(which(x<0)))
down.pace=down.pace[order(down.pace$diff,-down.pace$order1,-down.pace$order2,decreasing = T),]
heat.df=rbind(up.pace,down.pace)
heat.df[is.na(heat.df)] <- 0
write.csv(heat.df,"E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/all/Data_for_ADA_Heamap.csv",row.names = T)

heat.df <- read.csv("E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/all_sample/Data_for_ADA_Heamap.csv",row.names = 1)
max(heat.df)
min(heat.df)
p=pheatmap(heat.df,
           cluster_rows = F,
           cluster_cols =F,
           #scale =  "row", #"row""column"
           color =c(colorRampPalette(c("#2166AC","white"))(100),colorRampPalette(c("white","#E61C1D"))(100)) ,# "#2166AC" "#F7F7F7" "#B2182B"
           #color =c(colorRampPalette(c("grey95","#E6949A"))(500)) ,# "#2166AC" "#F7F7F7" "#B2182B"
           border_color = NA,
           #border_color = '#000000',
           fontsize = 14,
           fontsize_row = 10, #fontsize,
           fontsize_col = 10, #fontsize,
           fontsize_number = 0.8* fontsize,
           kmeans_k = NA,
           breaks = seq(-10 ,10,length.out =200),
           #treeheight_row=0,
           #cutree_cols = 2,
           #cutree_rows = 3,
           show_rownames = F,
           show_colnames = T,###列名是否保留
           legend = TRUE,#是不是加图例
           legend_breaks = NA,
           legend_labels = NA,
           drop_levels = TRUE,
           gaps_row = c(24,80),
           #gaps_col = seq(1,ncol(heat)),
           labels_row = NULL,
           labels_col = NULL,
           cellwidth = 10,
           cellheight = 1.5)#,#10 rescue
ggsave('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/all_sample/accelerators-decelerators_top20.pdf', plot = p, width =18, height =8)

heat <- dcast(df, sample+Age+gender ~ subtype, value.var = "group")
heat <- merge(heat, num,by.x = "sample", by.y = "Var1",all.x = T)
heat <- heat[,c('sample','Age','gender',"Group")]
rownames(heat) <- heat$sample
heat <- heat[rownames(heat.df),]
write.csv(heat,'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/all/accelerators-decelerators_top20.csv')

age <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/all/accelerators-decelerators_top20.csv',row.names = 1)
age <- age[,c('Age',"gender")]
p1=pheatmap(age,
           cluster_rows = F,
           cluster_cols =F,
           color =c(colorRampPalette(c("#13424c","#218490","#ae7c54"))(200)) ,# "#2166AC" "#F7F7F7" "#B2182B"
           border_color = NA,
           #border_color = '#000000',
           fontsize = 14,
           fontsize_row = 10, #fontsize,
           fontsize_col = 10, #fontsize,
           fontsize_number = 0.8* fontsize,
           kmeans_k = NA,
           show_rownames = F,
           show_colnames = T,###列名是否保留
           legend = TRUE,#是不是加图例
           legend_breaks = NA,
           legend_labels = NA,
           drop_levels = TRUE,
           breaks = seq(20 ,82,length.out =200),
           gaps_row = c(24,80),
           gaps_col = seq(1,ncol(age)),
           labels_row = NULL,
           labels_col = NULL,
           cellwidth = 10,
           cellheight = 1.5)#,#10 rescue
p2=pheatmap(age,
            cluster_rows = F,
            cluster_cols =F,
            color =c(colorRampPalette(c("#13424c","#218490","#ae7c54"))(200)) ,# "#2166AC" "#F7F7F7" "#B2182B"
            border_color = NA,
            #border_color = '#000000',
            fontsize = 14,
            fontsize_row = 10, #fontsize,
            fontsize_col = 10, #fontsize,
            fontsize_number = 0.8* fontsize,
            kmeans_k = NA,
            show_rownames = F,
            show_colnames = T,###列名是否保留
            legend = TRUE,#是不是加图例
            legend_breaks = NA,
            legend_labels = NA,
            drop_levels = TRUE,
            breaks = seq(0 ,1,length.out =200),
            gaps_row = c(24,80),
            gaps_col = seq(1,ncol(age)),
            labels_row = NULL,
            labels_col = NULL,
            cellwidth = 10,
            cellheight = 1.5)#,#10 rescue
p <- plot_grid(p1$gtable, p2$gtable,ncol = 2,align = "v")
ggsave('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/all/accelerators-decelerators_top20_age_gender.pdf', plot = p, width =18, height =8)


#####midutu#####
data <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/Age_pace_all.csv')
data <- subset(data, subtype == 'all_sample')

# 计算20%和80%分位数
q20 <- quantile(data$Age_pace, 0.2)
q80 <- quantile(data$Age_pace, 0.8)

p <- ggplot(data, aes(x = Age_pace)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = c(q20, q80), 
             color = "red", 
             linetype = "dashed",
             linewidth = 1) +
  labs(title = "密度图与20%/80%分位数线",
       x = "值",
       y = "密度") +
  theme_minimal()
ggsave('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/02_pace/all_sample_midutu.pdf', plot = p, width =8, height =6)
