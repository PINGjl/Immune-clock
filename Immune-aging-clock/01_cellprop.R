#install.packages('glmnet')
library(glmnet)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(cowplot)
library(hstats)

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

valid.rows <- sampled_prop$donor
train.rows <- setdiff(inf$donor, sampled_prop$donor)

########module#####
inf <- read.csv("E:/01_progrom/06_qz_pbmc_new/00_inf/inf_merge_final_me.csv")
cellprop <- read.csv("E:/01_progrom/06_qz_pbmc_new/03_celltype/final_celltype_prop.csv",row.names = 1)
cellprop <- merge(inf, cellprop, by="donor")

sign <- read.csv("E:/01_progrom/06_qz_pbmc_new/03_celltype/final_celltype_pcor.csv")
sign <- subset(sign, p.value < 0.05)
cellprop <- subset(cellprop, subtype %in% sign$subtype)

wtprop <- subset(cellprop, group  == 'control')
hrtprop <- subset(cellprop, group  == 'Used')
hrtprop = dcast(hrtprop, donor+gender ~ subtype, value.var = 'prop')
rownames(hrtprop) <- hrtprop$donor
hrtprop <- hrtprop[,-1]
hrtprop$gender <- as.numeric(as.factor(hrtprop$gender)) -1 ##female=0,male=1
inf <- inf[,c('donor', 'Age', 'gender')]

result <- 'E:/01_progrom/06_qz_pbmc_new/04_clock_gender/01_prop/result/'
mae.df=c()
for (alp in seq(0.01,0.99,0.01)){
set.seed(2024)
dir.create(paste0(result,alp))
pd <- unique(wtprop[,c('donor','Age', 'gender')])
pd$dataset <- pd$donor
pd[pd$dataset %in% train.rows,]$dataset <- 'Train'
pd[pd$dataset %in% valid.rows,]$dataset <- 'Valid'

shannodata= dcast(cellprop, donor+gender ~ subtype, value.var = 'prop')
rownames(shannodata) <- shannodata$donor
shannodata <- shannodata[,-1]
shannodata$gender <- as.numeric(as.factor(shannodata$gender)) -1 ##female=0,male=1
x=shannodata[train.rows,]

y = subset(pd, donor %in% train.rows)
rownames(y) <- y$donor
y <- y[rownames(x),]
y = as.numeric(y$Age)

fitcv <- cv.glmnet(as.matrix(x),y,alpha= alp)
model = glmnet(as.matrix(x),y, family="gaussian", alpha=alp, lambda = fitcv$lambda.min)
save(model, file= paste0(result,alp,"/model.RData"))

coef <- as.matrix(coef(model,s = "lambda.min"))
write.csv(coef,paste0(result,alp,'/prop_coef.csv'))

s <- perm_importance(model, y = y, X = as.matrix(x))
plot(s[1:10])+theme_bw()
ggsave(paste0(result,alp,'/perm_importance.pdf'),width =2.5,height = 2.5)
write.csv(s$M,paste0(result,alp,'/importance.csv'))
write.csv(s$SE,paste0(result,alp,'/importance_SE.csv'))

pred.train = as.numeric(predict(model, newx = as.matrix(x), s = "lambda.min"))
pred.val = as.numeric(predict(model, newx = as.matrix(shannodata[valid.rows,]), s = "lambda.min"))
train.res = data.frame(predAge = pred.train,row.names=train.rows)
val.res = data.frame(predAge = pred.val,row.names=valid.rows)
group=data.frame(donor=c(train.rows,valid.rows),pred_age = as.numeric(c(pred.train, pred.val)))
group=merge(group,pd,by='donor')
group=group[,c('donor','pred_age','Age', 'gender','dataset')]
write.csv(group,paste0(result,alp,'/pred_age.csv'))

pred.hrt = as.numeric(predict(model, newx = as.matrix(hrtprop), s = "lambda.min"))
hrt.res = data.frame(predAge = pred.hrt,donor=rownames(hrtprop))
pred.hrt=merge(inf, hrt.res, by = 'donor')
write.csv(group,paste0(result,alp,'/pred_age_hrt.csv'))

tmp = group
cor.train = rcorr(tmp[tmp$dataset=='Train','Age'], tmp[tmp$dataset=='Train','pred_age'], type = 'pearson') ###rcorr
cor1 = cor.train$r[2]
P1 = cor.train$P[2] 
cor.val = rcorr(tmp[tmp$dataset=='Valid','Age'], tmp[tmp$dataset=='Valid','pred_age'], type = 'pearson')
cor2 = cor.val$r[2]
P2 = cor.val$P[2]

mae1 = mean(abs(tmp[tmp$dataset=='Train','pred_age']-tmp[tmp$dataset=='Train','Age']))
mae2 = mean(abs(tmp[tmp$dataset=='Valid','pred_age']-tmp[tmp$dataset=='Valid','Age']))
mae.df=rbind(mae.df,data.frame(alpha=alp,MAE=mae2,R=cor2, P=P2))

train.plot = group[group$dataset=='Train',]
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
val.plot = group[group$dataset=='Valid',]
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

p <- plot_grid(p1, p2, ncol = 2,align = "v")
ggsave(paste0(result,alp,'/pred_age_module.pdf'), plot = p, width = 8, height =4)#7
}
mae.df <- mae.df[order(mae.df$MAE),]
write.csv(mae.df,paste0(result,'mae_result.csv'))
