library(immunarch)

data <- read.csv('E:/01_progrom/06_qz_pbmc_new/08_TBcr/tcr_meta.csv')
data <- subset(data, group == 'control')

freq <- as.data.frame(table(data$donor,data$CTaa))
freq <- subset(freq, Freq > 0)
data <- merge(data, freq, by.x = c('donor','CTaa'), by.y = c('Var1','Var2'))

num <- as.data.frame(table(data$donor))
colnames(num) <- c('donor','total')

data <- merge(data, num)
data$prop <- data$Freq / data$total

meta <- unique(data[,c('donor','Age')])
sample <- meta$donor
meta$donors <- c(paste0('D_',c(1:length(sample))))
samples <- meta$donors
for (i in (1:length(sample))){
  tmp <- subset(data, donor == sample[i])
  tmp <- unique(tmp[,c('CTaa','Freq','prop')])
  colnames(tmp) <- c('CDR3.aa','Clones','Proportion')
  assign(samples[i], tmp)
}

data <- list(D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8,D_9,D_10,D_11,D_12,D_13,D_14,D_15,D_16,D_17,D_18,D_19,D_20,D_21,D_22,D_23,D_24,D_25,D_26,D_27,D_28,D_29,D_30,D_31,D_32,D_33,D_34,D_35,D_36,D_37,D_38,D_39,D_40,D_41,D_42,D_43,
             D_44,D_45,D_46,D_47,D_48,D_49,D_50,D_51,D_52,D_53,D_54,D_55,D_56,D_57,D_58,D_59,D_60,D_61,D_62,D_63,D_64,D_65,D_66,D_67,D_68,D_69,D_70,D_71,D_72,D_73,D_74,D_75,D_76,D_77,D_78,D_79,D_80,D_81,D_82,D_83,D_84,
             D_85,D_86,D_87,D_88,D_89,D_90,D_91,D_92,D_93,D_94,D_95,D_96,D_97,D_98,D_99,D_100,D_101,D_102,D_103,D_104,D_105,D_106,D_107,D_108,D_109,D_110,D_111,D_112,D_113,D_114,D_115,D_116,D_117,D_118,D_119,D_120,D_121,
             D_122,D_123,D_124,D_125,D_126,D_127,D_128,D_129,D_130,D_131,D_132,D_133,D_134,D_135,D_136,D_137,D_138,D_139,D_140,D_141,D_142,D_143,D_144,D_145,D_146,D_147,D_148,D_149,D_150,D_151,D_152,D_153,D_154,D_155,
             D_156,D_157,D_158,D_159,D_160,D_161,D_162,D_163,D_164,D_165,D_166,D_167,D_168,D_169,D_170,D_171,D_172,D_173,D_174,D_175,D_176,D_177,D_178,D_179,D_180,D_181,D_182,D_183,D_184,D_185,D_186,D_187,D_188,D_189,
             D_190,D_191,D_192,D_193)
names(data) <- samples

immdata <- list(meta, data)
names(immdata) <- c('meta','data')
write.csv(meta,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/meta_inf.csv',row.names = F)

# Compute statistics and visualise them
# Chao1 diversity measure
div_chao <- repDiversity(immdata$data, "chao1")
write.csv(div_chao,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/Chao1.csv')
# Hill numbers
div_hill <- repDiversity(immdata$data, "hill")
write.csv(div_hill,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/Hill_numbers.csv')
# Gini-Simpson index
div_ginisimp <- repDiversity(immdata$data, "gini.simp")
write.csv(div_ginisimp,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/Gini_Simpson.csv')
# Gini coefficient
div_gini <- repDiversity(immdata$data, "gini")
write.csv(div_gini,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/Gini.csv')
# Inverse Simpson index
div_Inverse <- repDiversity(immdata$data, "inv.simp")
write.csv(div_Inverse,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/Inverse_Simpson.csv')
# D50
div_d50 <- repDiversity(immdata$data, "d50")
write.csv(div_d50,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/D50.csv')
# Ecological diversity measure
div_div <- repDiversity(immdata$data, "div")
write.csv(div_div,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/Ecological_diversity_measure.csv')
## Rarefaction 文件太大，没用到
#imm_raref <- repDiversity(immdata$data, "raref", .verbose = F)
#write.csv(imm_raref,'E:/01_progrom/06_qz_pbmc_new/08_TBcr/02_immunarch/Rarefaction.csv')
