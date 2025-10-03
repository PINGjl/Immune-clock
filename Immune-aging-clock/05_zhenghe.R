####mae/plot#####
celltypes <- c("Tn_CD4", "Tcm_CD4","Th2_CD4","Th17_CD4","Treg",
               "Tn_CD8","Tem","TCR","CTL",
               "Tex_CD8",'proT','NK1','NK2','Bm',
               'pDC','CD4T','CD8T',"cMC","Pla","HSPC")
celltypes <- c("cDC") #97
celltypes <- c('Tcm_CD8') #99
celltypes <- c("Th1_CD4") #-103
celltypes <- c('nMC') #3
celltypes <- c("Bn") #0-101
celltypes <- c("Tex_CD4") #1-  alpha = 0.02 ——重跑重要性
#celltypes <- c('cMC') #qu diao 	male-pbmc-15 ——重跑年龄预测


mae.df <- c()
df <- c()
for (cell in celltypes) {
    result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/02_gene/',cell,'/')
  
    mae=read.csv(paste0(result,'mae_result.csv'))
    mae=mae[order(mae$MAE),]
    j=mae$alpha[1]
    
    pred <- read.csv(paste0(result,j,'/pred_age.csv'))
    over <- subset(pred, pred_age < 10 | pred_age > 95)
    # 如果 over 有数据，才合并
    if (nrow(over) > 0) {
      df <- rbind(df, cbind(over, subtype = cell))
    }
    
    cor.train = rcorr(pred[pred$dataset=='Train','Age'], pred[pred$dataset=='Train','pred_age'], type = 'pearson') ###rcorr
    cor1 = cor.train$r[2]
    P1 = cor.train$P[2] 
    cor.val = rcorr(pred[pred$dataset=='Valid','Age'], pred[pred$dataset=='Valid','pred_age'], type = 'pearson')
    cor2 = cor.val$r[2]
    P2 = cor.val$P[2]
    
    mae1 = mean(abs(pred[pred$dataset=='Train','pred_age']-pred[pred$dataset=='Train','Age']))
    mae2 = mean(abs(pred[pred$dataset=='Valid','pred_age']-pred[pred$dataset=='Valid','Age']))
    mae.df=rbind(mae.df,data.frame(subtype=cell,alpha=j,MAE=mae2,R=cor2, P=P2))
    
    train.plot = pred[pred$dataset=='Train',]
    p1 = ggplot() +
      geom_point(data = train.plot, aes(Age, pred_age), size = 2) +
      geom_abline(data = NULL,linetype=4,color="grey") +
      geom_smooth(data = train.plot, aes(Age, pred_age), method = lm, se = F, color = '#3d60b7') +
      geom_text(aes(x = 15, y = 70), hjust = 0,
                label = paste0('R = ', round(cor1, 2), '\nP = ', P1, '\nMAE = ', round(mae1, 2))) +
      scale_y_continuous(limits = c(10, 95)) +
      scale_x_continuous(limits = c(10, 95)) +
      labs(y = paste0('Train')) +
      theme_base() +
      theme(legend.position = 'none')
    val.plot = pred[pred$dataset=='Valid',]
    p2 = ggplot() +
      geom_point(data = val.plot, aes(Age, pred_age), size = 2) +
      geom_abline(data = NULL,linetype=4,color="grey") +
      geom_smooth(data = val.plot, aes(Age, pred_age), method = lm, se = F, color = '#3d60b7') +
      geom_text(aes(x = 15, y = 70), hjust = 0,
                label = paste0('R = ', round(cor2, 2), '\nP = ', P2, '\nMAE = ', round(mae2, 2))) +
      scale_y_continuous(limits = c(1, 95)) +
      scale_x_continuous(limits = c(1, 95)) +
      labs(y = paste0('Valid')) +
      theme_base() +
      theme(legend.position = 'none')
    
    p <- plot_grid(p1, p2,ncol = 2,align = "v")
    ggsave(paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/01_plot/',cell,'_',j,'_pred_age_module_2.pdf'), plot = p, width = 8, height =4)#7
}

celltypes <- c('tc','bc','nk','mc','dc')#'all','all_sample'
for (cell in celltypes) {
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/04_combind/',cell,'/')
  
  mae=read.csv(paste0(result,'mae_result.csv'))
  mae=mae[order(mae$MAE),]
  j=mae$alpha[1]
  
  pred <- read.csv(paste0(result,j,'/pred_age.csv'))
  over <- subset(pred, pred_age < 10 | pred_age > 95)
  # 如果 over 有数据，才合并
  if (nrow(over) > 0) {
    df <- rbind(df, cbind(over, subtype = cell))
  }
  
  cor.train = rcorr(pred[pred$dataset=='Train','Age'], pred[pred$dataset=='Train','pred_age'], type = 'pearson') ###rcorr
  cor1 = cor.train$r[2]
  P1 = cor.train$P[2] 
  cor.val = rcorr(pred[pred$dataset=='Valid','Age'], pred[pred$dataset=='Valid','pred_age'], type = 'pearson')
  cor2 = cor.val$r[2]
  P2 = cor.val$P[2]
  
  mae1 = mean(abs(pred[pred$dataset=='Train','pred_age']-pred[pred$dataset=='Train','Age']))
  mae2 = mean(abs(pred[pred$dataset=='Valid','pred_age']-pred[pred$dataset=='Valid','Age']))
  mae.df=rbind(mae.df,data.frame(subtype=cell,alpha=j,MAE=mae2,R=cor2, P=P2))
  
  train.plot = pred[pred$dataset=='Train',]
  p1 = ggplot() +
    geom_point(data = train.plot, aes(Age, pred_age), size = 2) +
    geom_abline(data = NULL,linetype=4,color="grey") +
    geom_smooth(data = train.plot, aes(Age, pred_age), method = lm, se = F, color = '#3d60b7') +
    geom_text(aes(x = 15, y = 70), hjust = 0,
              label = paste0('R = ', round(cor1, 2), '\nP = ', P1, '\nMAE = ', round(mae1, 2))) +
    scale_y_continuous(limits = c(10, 95)) +
    scale_x_continuous(limits = c(10, 95)) +
    labs(y = paste0('Train')) +
    theme_base() +
    theme(legend.position = 'none')
  val.plot = pred[pred$dataset=='Valid',]
  p2 = ggplot() +
    geom_point(data = val.plot, aes(Age, pred_age), size = 2) +
    geom_abline(data = NULL,linetype=4,color="grey") +
    geom_smooth(data = val.plot, aes(Age, pred_age), method = lm, se = F, color = '#3d60b7') +
    geom_text(aes(x = 15, y = 70), hjust = 0,
              label = paste0('R = ', round(cor2, 2), '\nP = ', P2, '\nMAE = ', round(mae2, 2))) +
    scale_y_continuous(limits = c(10, 95)) +
    scale_x_continuous(limits = c(10, 95)) +
    labs(y = paste0('Valid')) +
    theme_base() +
    theme(legend.position = 'none')
  
  p <- plot_grid(p1, p2,ncol = 2,align = "v")
  ggsave(paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/01_plot/',cell,'_',j,'_pred_age_module_2.pdf'), plot = p, width = 8, height =4)#7
}
for (cell in c('result')) {
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/01_prop/',cell,'/') #03_tcr
  
  mae=read.csv(paste0(result,'mae_result.csv'))
  mae=mae[order(mae$MAE),]
  j=mae$alpha[1]
  
  pred <- read.csv(paste0(result,j,'/pred_age.csv'))
  over <- subset(pred, pred_age < 10 | pred_age > 95)
  # 如果 over 有数据，才合并
  if (nrow(over) > 0) {
    df <- rbind(df, cbind(over, subtype = cell))
  }
  
  cor.train = rcorr(pred[pred$dataset=='Train','Age'], pred[pred$dataset=='Train','pred_age'], type = 'pearson') ###rcorr
  cor1 = cor.train$r[2]
  P1 = cor.train$P[2] 
  cor.val = rcorr(pred[pred$dataset=='Valid','Age'], pred[pred$dataset=='Valid','pred_age'], type = 'pearson')
  cor2 = cor.val$r[2]
  P2 = cor.val$P[2]
  
  mae1 = mean(abs(pred[pred$dataset=='Train','pred_age']-pred[pred$dataset=='Train','Age']))
  mae2 = mean(abs(pred[pred$dataset=='Valid','pred_age']-pred[pred$dataset=='Valid','Age']))
  mae.df=rbind(mae.df,data.frame(subtype=cell,alpha=j,MAE=mae2,R=cor2, P=P2))
  
  train.plot = pred[pred$dataset=='Train',]
  p1 = ggplot() +
    geom_point(data = train.plot, aes(Age, pred_age), size = 2) +
    geom_abline(data = NULL,linetype=4,color="grey") +
    geom_smooth(data = train.plot, aes(Age, pred_age), method = lm, se = F, color = '#3d60b7') +
    geom_text(aes(x = 15, y = 70), hjust = 0,
              label = paste0('R = ', round(cor1, 2), '\nP = ', P1, '\nMAE = ', round(mae1, 2))) +
    scale_y_continuous(limits = c(10, 95)) +
    scale_x_continuous(limits = c(10, 95)) +
    labs(y = paste0('Train')) +
    theme_base() +
    theme(legend.position = 'none')
  val.plot = pred[pred$dataset=='Valid',]
  p2 = ggplot() +
    geom_point(data = val.plot, aes(Age, pred_age), size = 2) +
    geom_abline(data = NULL,linetype=4,color="grey") +
    geom_smooth(data = val.plot, aes(Age, pred_age), method = lm, se = F, color = '#3d60b7') +
    geom_text(aes(x = 15, y = 70), hjust = 0,
              label = paste0('R = ', round(cor2, 2), '\nP = ', P2, '\nMAE = ', round(mae2, 2))) +
    scale_y_continuous(limits = c(10, 95)) +
    scale_x_continuous(limits = c(10, 95)) +
    labs(y = paste0('Valid')) +
    theme_base() +
    theme(legend.position = 'none')
  
  p <- plot_grid(p1, p2,ncol = 2,align = "v")
  ggsave(paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/01_plot/',cell,'_',j,'_pred_age_module_2.pdf'), plot = p, width = 8, height =4)#7
}

mae.df=mae.df[order(mae.df$MAE),]
write.csv(mae.df,'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/mae_all.csv')
write.csv(df,'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/over_predAge.csv')

mae.df <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/mae_all.csv')
mae.df <- subset(mae.df, !subtype %in% c('CD4T','CD8T', 'all'))
p <- ggplot(mae.df, aes(R, (1/MAE))) +
  geom_point(aes(color=Group),size =3)+ #size =R, 
  theme_base()
ggsave('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/mae_all.pdf', plot = p, width = 6.1, height =4.5)

###predAge####
mae.df <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/mae_all.csv')
df <- c()

gene <- subset(mae.df, Group == 'tAge')
for (cell in unique(gene$subtype)) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/02_gene/',cell,'/',tmp$alpha,'/')
  
  pred <- read.csv(paste0(result,'pred_age.csv'),row.names = 1)
  pred <- subset(pred, dataset == 'Valid')
  pred$subtype <- cell
  
  df=rbind(df,pred)
}

gene <- subset(mae.df, Group %in% c('iAge','ptAge'))
for (cell in unique(gene$subtype)) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/04_combind/',cell,'/',tmp$alpha,'/')
  
  pred <- read.csv(paste0(result,'pred_age.csv'),row.names = 1)
  pred <- subset(pred, dataset == 'Valid')
  pred$subtype <- cell
  
  df=rbind(df,pred)
}

gene <- subset(mae.df, Group == 'pAge')
for (cell in c('result')) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/01_prop/',cell,'/',tmp$alpha,'/')
  
  pred <- read.csv(paste0(result,'pred_age.csv'),row.names = 1)
  pred <- subset(pred, dataset == 'Valid')
  pred$subtype <- 'pAge'
  
  df=rbind(df,pred)
}

gene <- subset(mae.df, Group == 'TCRAge')
for (cell in c('result')) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/03_tcr/',cell,'/',tmp$alpha,'/')
  
  pred <- read.csv(paste0(result,'pred_age.csv'),row.names = 1)
  pred <- subset(pred, dataset == 'Valid')
  pred$subtype <- 'TCRAge'
  
  df=rbind(df,pred)
}

write.csv(df,'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/predAge_all.csv', row.names = F)


###coefficient importance####
mae.df <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/mae_all.csv')
df <- c()

gene <- subset(mae.df, Group == 'tAge')
for (cell in unique(gene$subtype)) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/02_gene/',cell,'/',tmp$alpha,'/')
  
  coef <- read.csv(paste0(result,'prop_coef.csv'))
  coef <- subset(coef, !s1 == 0)
  
  impor <- read.csv(paste0(result,'importance.csv'))
  impor <- subset(impor, !s0 == 0)
  se <- read.csv(paste0(result,'importance_SE.csv'))
  se <- subset(se, !s0 == 0)
  
  final <- merge(coef, impor, by = 'X', all=T)
  final <- merge(final, se, by = 'X', all=T)
  colnames(final) <- c('Feature','Coef','Importance','Importance_SE')
  final$subtype <- cell
  
  df=rbind(df,final)
}

gene <- subset(mae.df, Group %in% c('iAge','ptAge'))
for (cell in unique(gene$subtype)) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/04_combind/',cell,'/',tmp$alpha,'/')
  
  coef <- read.csv(paste0(result,'prop_coef.csv'))
  coef <- subset(coef, !s1 == 0)
  
  impor <- read.csv(paste0(result,'importance.csv'))
  impor <- subset(impor, !s0 == 0)
  se <- read.csv(paste0(result,'importance_SE.csv'))
  se <- subset(se, !s0 == 0)
  
  final <- merge(coef, impor, by = 'X', all=T)
  final <- merge(final, se, by = 'X', all=T)
  colnames(final) <- c('Feature','Coef','Importance','Importance_SE')
  final$subtype <- cell
  
  df=rbind(df,final)
}

gene <- subset(mae.df, Group == 'pAge')
for (cell in c('result')) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/01_prop/',cell,'/',tmp$alpha,'/')
  
  coef <- read.csv(paste0(result,'prop_coef.csv'))
  coef <- subset(coef, !s1 == 0)
  
  impor <- read.csv(paste0(result,'importance.csv'))
  impor <- subset(impor, !s0 == 0)
  se <- read.csv(paste0(result,'importance_SE.csv'))
  se <- subset(se, !s0 == 0)
  
  final <- merge(coef, impor, by = 'X', all=T)
  final <- merge(final, se, by = 'X', all=T)
  colnames(final) <- c('Feature','Coef','Importance','Importance_SE')
  final$subtype <- 'pAge'
  
  df=rbind(df,final)
}

gene <- subset(mae.df, Group == 'TCRAge')
for (cell in c('result')) {
  tmp <- subset(gene, subtype == cell)
  result <- paste0('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/03_tcr/',cell,'/',tmp$alpha,'/')
  
  coef <- read.csv(paste0(result,'prop_coef.csv'))
  coef <- subset(coef, !s1 == 0)
  
  impor <- read.csv(paste0(result,'importance.csv'))
  impor <- subset(impor, !s0 == 0)
  se <- read.csv(paste0(result,'importance_SE.csv'))
  se <- subset(se, !s0 == 0)
  
  final <- merge(coef, impor, by = 'X', all=T)
  final <- merge(final, se, by = 'X', all=T)
  colnames(final) <- c('Feature','Coef','Importance','Importance_SE')
  final$subtype <- 'TCRAge'
  
  df=rbind(df,final)
}

write.csv(df,'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/importance_all.csv', row.names = F)

######importance plot#####
df <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/importance_all.csv')

result <- subset(df, subtype == "all_sample") #TCRAge
result <- result[-1,]
result[result$Feature %in% c('Hill_numbers_Q5','Hill_numbers_Q6'),]$subtype <- 'TCR'
result[result$Feature %in% c('Tn_CD8'),]$subtype <- 'prop'

result <- result[order(result$Importance),]
result$Feature <- factor(result$Feature, 
                      levels=result$Feature, ordered=TRUE) 
value <- max(result$Importance)
result$Importance <- result$Importance/ value

p1 <- ggplot(result, aes(x = Importance, y = Feature, fill = subtype)) +
  geom_bar(stat = "identity", position = "dodge") 

result <- result[order(result$Importance, decreasing = T),]
value <- max(result$Importance)
result$Importance <- result$Importance/ value
result$Importance_SE <- result$Importance_SE / value
importance <- result#[c(1:10),]
###具体指标
p2 <- ggplot(importance, aes(x=Importance, y=Feature, colour=subtype)) + 
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(xmin=Importance-Importance_SE, xmax=Importance+Importance_SE), width=.3) +
  #geom_line() +
  #geom_point()+
  #scale_color_manual(values=CT.col) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'),
        axis.text.y = element_text(color='black'),
        axis.text.x = element_text(color = 'black'))#, angle = 90, hjust=1, vjust=1))
p <- plot_grid(p1, p2, ncol = 2,align = "v")
ggsave("E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/importance_pAge.pdf", plot = p2, width = 7, height =5)#7

######gene plot#####
celltypes <- c("Tn_CD4", "Tcm_CD4","Th1_CD4","Th2_CD4","Th17_CD4","Treg","Tex_CD4",
               "Tn_CD8",'Tcm_CD8',"Tem","TCR","CTL","Tex_CD8",'proT',
               'NK1','NK2','Bn','Bm','Pla','cMC','nMC','cDC','pDC','HSPC')
df <- read.csv('E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/importance_all.csv')
df <- subset(df, subtype %in% celltypes)
df <- subset(df, ! Importance == '')

df <- df[order(df$Importance),]
top_100 <- df %>%
  group_by(subtype) %>%
  slice_max(order_by = Importance,n = 100,with_ties = TRUE)
meta <- dcast(top_100,Feature ~ subtype, value.var = 'Feature')
meta <- meta[,-1]
write.csv(meta,'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/transcript_importance_top_100_meta.csv', row.names = F)

num <- as.data.frame(table(top_100$Feature))
num$type <- as.factor(ifelse(num$Freq > 1 ,'Common','Specific'))
top_100 <- merge(top_100, num, by.x = 'Feature', by.y = 'Var1')
num <- as.data.frame(table(top_100$subtype, top_100$type))
p <- ggplot(num, aes(x = "", y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", width = 1) + # 使用geom_bar并设置宽度为1
  coord_polar(theta = "y") +
  facet_grid(.~ Var1) +theme_void()
pdf(paste("E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/transcript_importance_top_100_meta_pie.pdf",sep=""),width=10, height=2)
print(p)
dev.off()

df <- subset(df, ! Coef == '')
num <- as.data.frame(table(df$subtype))
num$Freq <- num$Freq -1
num <- num[order(num$Freq),]

num$Var1 <- factor(num$Var1, 
                      levels=num$Var1, ordered=TRUE) 
p <- ggplot(num, aes(x = Freq, y = Var1)) +
  geom_bar(stat = "identity", width = 0.7) #+ # 使用geom_bar并设置宽度为1
  #coord_polar(theta = "y") +
  #facet_grid(.~ Var1) +theme_void()
pdf(paste("E:/01_progrom/06_qz_pbmc_new/04_clock_gender/05_zhenghe/transcript_coef_num.pdf",sep=""),width=10, height=2)
print(p)
dev.off()
