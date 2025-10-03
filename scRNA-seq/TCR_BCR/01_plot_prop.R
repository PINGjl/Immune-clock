#####¼ì²â±ÈÀı#####
#TC 898015; BC 89908
tcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/tcr_meta.csv')
bcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/bcr_meta.csv')

######Diversity####
data <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/tcr_meta.csv')
data <- subset(data, group == 'control')

data2 <- unique(data[,c('donor','Age','CTaa')])
num <- as.data.frame(table(data2$donor,data2$Age)) #,data2$celltype
num <- subset(num, Freq > 0)
colnames(num) <- c('donor','Age','Diversity') #,'celltype'

cellnum <- as.data.frame(table(data$donor))#,data$celltype
colnames(cellnum) <- c('donor','Total') #,'celltype'

num <- merge(num, cellnum, by = 'donor')# ,'celltype','Var2'
num$Diversity_prop <- (num$Diversity / num$Total )*100
num$Age <- as.numeric(as.character(num$Age))
write.csv(num, "E:/01_progrom/06_qz_pbmc_new/08_TBcr/00_clock_inf/Diversity_donor_T.csv",row.names = F)

p1 = ggplot(num, aes(Age, Diversity_prop))+
  geom_point() +
  geom_smooth(method = "lm")+ 
  stat_cor(data=num, label.x.npc = "left",label.y.npc = "top",method = "pearson")+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))
p <- plot_grid(p1, p2, ncol = 2,align = "v")
pdf("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/TB_Diversity_pre_100.pdf",width=6, height=3)
print(p)
dev.off()

###celltype
tcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/tcr_meta.csv')
bcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/bcr_meta.csv')
data <- rbind(tcr,bcr)
data <- subset(data, group == 'control')

data2 <- unique(data[,c('donor','Age',"subtype",'CTaa')])
num <- as.data.frame(table(data2$donor,data2$Age,data2$subtype)) 
num <- subset(num, Freq > 0)
colnames(num) <- c('donor','Age','subtype','Diversity') 

cellnum <- as.data.frame(table(data$donor,data$subtype))
colnames(cellnum) <- c('donor','subtype','Total') 

num <- merge(num, cellnum, by = c('donor','subtype'))# ,'celltype','Var2'
num$Diversity_prop <- (num$Diversity / num$Total )*100
num$Age <- as.numeric(as.character(num$Age))

df <- data.frame()
for (i in unique(num$subtype)){
  obj <- subset(num, subtype == i)
  celltype <- i
  cor = rcorr(obj$Age,obj$Diversity_prop)
  coef = cor$r[2]
  P = cor$P[2] 
  result <- data.frame(celltype,coef,P)
  df <- rbind(df,result)
}
df$logP <- -(log10(df$P))
df <- df[order(df$logP),]#, decreasing = T
df$celltype <- factor(df$celltype, 
                      levels=c(df$celltype), ordered=TRUE) 

num$subtype <- factor(num$subtype, 
                      levels=rev(df$celltype), ordered=TRUE)
p = ggplot(num, aes(Age, Diversity_prop))+
  facet_wrap(~subtype, scales = "free",ncol = 6)+
  geom_point() +
  #scale_color_manual(values=cols) +
  geom_smooth(method = "lm")+ stat_cor(data=num, label.x.npc = "left",
                                       label.y.npc = "top",method = "pearson")+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))
pdf(paste("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/TB_celltype_Diversity_pre_100.pdf",sep=""),width=13.5, height=7.5)
print(p)
dev.off()

df <- subset(df, (coef < 0 & P < 0.05))
p = ggplot(df, aes(celltype, coef))+
  geom_point(size = 2)+ylim(-0.3,0)
pdf(paste("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/TB_celltype_Diversity_pre_100-2.pdf",sep=""),width=3, height=3)
print(p)
dev.off()

num <- reshape2::dcast(num, donor ~ subtype, value.var = "Diversity_prop")
write.csv(num, "E:/01_progrom/06_qz_pbmc_new/08_TBcr/00_clock_inf/Diversity_donor_T_subtype.csv",row.names = F)

######cloneType analysis#####
tcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/tcr_meta.csv')
tcr$anno <- 'TC'
bcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/bcr_meta.csv')
bcr$anno <- 'BC'
data <- rbind(tcr,bcr)
data <- subset(data, group == 'control')

freq <- as.data.frame(table(data$donor,data$anno,data$CTaa))
freq <- subset(freq, Freq > 0)

data <- merge(data, freq, by.x = c('donor','anno','CTaa'), by.y = c('Var1','Var2','Var3'))
data$cloneType <- as.factor(ifelse(data$Freq == 1,'Single','Multiple'))
#data <- subset(data, anno == 'TC')
#write.csv(data, 'E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/cloneType_donor_T-Multiple_Single.csv',row.names = F)

num <- as.data.frame(table(data$donor,data$anno,data$cloneType))
num <- subset(num, Freq > 0)
colnames(num) <- c('donor','anno','cloneType','number')
num <- num %>% group_by(donor, anno) %>% mutate(total=sum(number))
num$prop <- num$number / num$total

inf <- unique(data[,c('donor','Age')])
num <- merge(inf, num, by="donor")

num$cloneType <- factor(num$cloneType, 
                       levels=c("Single",'Multiple'), ordered=TRUE)

cols <- c('#004597','#881600')
p = ggplot(num, aes(Age, prop, colour = cloneType))+
  facet_wrap(~ anno, scales = "free",ncol = 2)+
  geom_point() +
  scale_color_manual(values=cols) +
  geom_smooth(method = "lm")+ stat_cor(data=num, label.x.npc = "left",
                                       label.y.npc = "top",method = "pearson")+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))
pdf("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/cloneType_donor_T-B_prop.pdf",width=6.9, height=3.2)
print(p)
dev.off()

num <- subset(num, anno == 'TC')
num <- dcast(num, donor~cloneType, value.var = 'prop')
write.csv(num, "E:/01_progrom/06_qz_pbmc_new/08_TBcr/00_clock_inf/cloneType_donor_T_prop.csv",row.names = F)

###celltype
data <- subset(data, anno == 'TC')
total <- as.data.frame(table(data$subtype))
colnames(total) <- c('subtype','total')
kuozeng <- subset(data, Freq > 1)
num <- as.data.frame(table(kuozeng$subtype))
num <- merge(total,num, by.x = 'subtype', by.y = 'Var1')
num$prop <- num$Freq / num$total
num$weight <- num$prop / sum(num$prop)
write.csv(num, "E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/Tcell_cloneType.csv",row.names = F)
num <- read.csv("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/Tcell_cloneType.csv")
num$h2 <- factor(num$h2, 
                        levels=num$h2, ordered=TRUE)

p <- ggplot(num, aes(h1,weight,fill = h2)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.8)+
  scale_fill_manual(values = num$color)
pdf("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/Tcell_cloneType.pdf", width =3, height =5)
print(p)
dev.off()

######max clones num####
tcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/tcr_meta.csv')
tcr$anno <- 'TC'
bcr <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/bcr_meta.csv')
bcr$anno <- 'BC'
data <- rbind(tcr,bcr)
data <- subset(data, group == 'control')

freq <- as.data.frame(table(data$donor,data$anno,data$CTaa))
freq <- subset(freq, Freq > 0)
data <- merge(data, freq, by.x = c('donor','anno','CTaa'), by.y = c('Var1','Var2','Var3'))
data$cloneType <- as.factor(ifelse(data$Freq == 1,'Single','Multiple'))

num <- as.data.frame(table(data$donor,data$anno))
colnames(num) <- c('donor','anno','total')
data <- merge(data, num, by = c('donor','anno'))
data$prop <- data$Freq / data$total

df <- unique(data[,c("donor", "anno","Age","Freq","total","prop")])
df <- df[order(df$prop,decreasing = T),]
df$ident <- paste0(df$donor,df$anno)
max <- df[!duplicated(df$ident),]

p = ggplot(max, aes(Age, prop))+
  facet_wrap(~anno, scales = "free",ncol = 2)+
  geom_point() +
  #scale_color_manual(values=cols) +
  geom_smooth(method = "lm")+ stat_cor(data=max, label.x.npc = "left",
                                       label.y.npc = "top",method = "pearson")+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))
pdf("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/TB_max_clones_prop.pdf",width=6, height=3.3)
print(p)
dev.off()

max <- subset(max, anno == 'TC')
colnames(max)[6] <- 'Max_clones_prop' 
max <- max[,c(1,6)]
write.csv(max, "E:/01_progrom/06_qz_pbmc_new/08_TBcr/00_clock_inf/cloneType_donor_T_Max_clones_prop.csv",row.names = F)

######Public or Private####
data <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/bcr_meta.csv')
data <- subset(data, group == 'control')

type <- unique(data[,c("donor","CTaa")])
num <- as.data.frame(table(type$CTaa))
colnames(num) <- c('CTaa','num')
num$type <- as.factor(ifelse(num$num <= 1,'Private','Public'))

data <- merge(data, num, by='CTaa')
num <- as.data.frame(table(data$donor,data$type))
colnames(num) <- c('donor','type','number')
num <- num %>% group_by(donor) %>% mutate(total=sum(number))
num$prop <- num$number / num$total

inf <- unique(data[,c('donor','Age')])
num <- merge(inf, num, by="donor")
#num <- dcast(num, donor~type, value.var = 'prop')
#write.csv(num, "E:/01_progrom/06_qz_pbmc_new/08_TBcr/00_clock_inf//cloneType_donor_T_Public_Private.csv",row.names = F)

p2 = ggplot(num, aes(Age, prop,color = type))+
  #facet_wrap(~type, scales = "free",ncol = 2)+
  geom_point() +
  scale_color_manual(values=c('#004C63','#EA5E33')) +
  geom_smooth(method = "lm")+ stat_cor(data=num, label.x.npc = "left",
                                       label.y.npc = "top",method = "pearson")+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'right',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))
p <- plot_grid(p1, p2, ncol = 2,align = "v")
pdf("E:/01_progrom/06_qz_pbmc_new/08_TBcr/01_prop/T-B_Public.pdf",width=8.2, height=3)
print(p)
dev.off()

