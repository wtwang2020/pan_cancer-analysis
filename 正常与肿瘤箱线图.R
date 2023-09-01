#绘制肿瘤与正常组织之间各基因的表达情况
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/图片区/神经肿瘤项目/单基因")
library(installr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
######################该基因集输入和输出名称,更改基因集时循环长度可能要改，长的话###################################
#读取数据
load("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/output_pancancer_xena/TCGA_GTEx_pancancer_mrna_pheno.rdata")
#genes <- c('ADCYAP1R1',	'ADIPOR1',	'ADIPOR2',	'AGTR1',	'AGTR2',	'APLNR',	'AVPR1A',	'AVPR1B',	'AVPR2',	'BDKRB1',	'BDKRB2',	'CALCR',	'CALCRL',	'CCKAR',	'CCKBR',	'CRHR1',	'CRHR2',	'EDNRA',	'EDNRB',	'GALR1',	'GALR2',	'GALR3',	'GCGR',	'GHR',	'GHRHR',	'GHSR',	'GIPR',	'GLP1R',	'GLP2R',	'GNRHR',	'GRPR',	'HCRTR1',	'HCRTR2',	'KISS1R',	'LEPR',	'MC3R',	'MC4R',	'MCHR1',	'MCHR2',	'MLNR',	'NMBR',	'NMUR1',	'NMUR2',	'NPBWR1',	'NPBWR2',	'NPFFR1',	'NPFFR2',	'NPR1',	'NPR2',	'NPR3',	'NPSR1',	'NPY1R',	'NPY2R',	'NPY4R',	'NPY5R',	'NTSR1',	'NTSR2',	'OPRD1',	'OPRK1',	'OPRL1',	'OPRM1',	'OXTR',	'PRLHR',	'PRLR',	'RAMP1',	'RAMP2',	'RAMP3',	'RXFP1',	'RXFP2',	'SCTR',	'SORT1',	'SSTR1',	'SSTR2',	'SSTR3',	'SSTR4',	'SSTR5',	'TACR1',	'TACR2',	'TACR3',	'TRHR',	'UTS2R',	'VIPR1',	'VIPR2')
genes <- c("RLN1", "GHRL", "CALCA", "PMCH", "PDYN", "UTS2B", "VIP", "GRP", "TAC1", "PENK",
          "SCG5", "CBLN4", "TAC4", "NPPC", "GNRH1", "EDN3", "EDN1", "CALCB", "NUCB2",
          "NPFF", "NMU", "RLN2", "KNG1", "CHGB", "CBLN2", "CARTPT", "POMC", "NPB", "UCN",
          "CORT", "AGT", "SCT", "HCRT", "AGRP", "VGF", "KISS1", "GAST", "UCN2", "GAL",
          "TAC3", "CRH", "SST", "ADCYAP1", "PNOC", "CBLN1", "CBLN3", "OXT", "NPW", "GNRH2",
          "ADM", "PYY", "NPPB", "TRH", "NPPA", "NMB", "NPY", "CHGA", "CCK", "GHRH", "APLN")
for (i in 1:83){
  print(genes[i])
gene <- genes[i]
plot_df <-tcga_gtex_mrna_pheno %>%
  select(1:4,all_of(gene)) %>%
  filter(sample_type %in% c("GTEx_normal", "TCGA_tumor"))
###############################################
p <- ggplot(plot_df, aes(project, get(gene))) +
  geom_boxplot(aes(fill = sample_type)) +
  theme(legend.position = "top") +
  labs(x = NULL, y = paste(gene, " expression")) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = plot_df$sample_type),label = "p.signif",color = "red")
#导出来图片
output_path <- "D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/图片区/神经肿瘤项目/单基因/配体"
output_filename <- paste(gene,".png")
# 设置图片的尺寸和分辨率
width <- 10  # 宽度（单位：英寸）
height <- 6  # 高度（单位：英寸）
dpi <- 1000   # 分辨率（dots per inch，每英寸的像素数）
p
# 使用 ggsave() 导出高清图片
ggsave(
  file = file.path(output_path, output_filename),
  plot = p,
  width = width,
  height = height,
  dpi = dpi
)
}
