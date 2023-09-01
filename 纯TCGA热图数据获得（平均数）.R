library(ggplot2)
library(tidyverse)
library(pheatmap)
#加载mRNA数据
load("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/TCGA纯癌症数据矩阵（基因已经选择）")
#####################输入基因集##########################
######################输入基因集（只需要改基因集合和输出文件）###################################
#gene <- c('ADCYAP1R1',	'ADIPOR1',	'ADIPOR2',	'AGTR1',	'AGTR2',	'APLNR',	'AVPR1A',	'AVPR1B',	'AVPR2',	'BDKRB1',	'BDKRB2',	'CALCR',	'CALCRL',	'CCKAR',	'CCKBR',	'CRHR1',	'CRHR2',	'EDNRA',	'EDNRB',	'GALR1',	'GALR2',	'GALR3',	'GCGR',	'GHR',	'GHRHR',	'GHSR',	'GIPR',	'GLP1R',	'GLP2R',	'GNRHR',	'GRPR',	'HCRTR1',	'HCRTR2',	'KISS1R',	'LEPR',	'MC3R',	'MC4R',	'MCHR1',	'MCHR2',	'MLNR',	'NMBR',	'NMUR1',	'NMUR2',	'NPBWR1',	'NPBWR2',	'NPFFR1',	'NPFFR2',	'NPR1',	'NPR2',	'NPR3',	'NPSR1',	'NPY1R',	'NPY2R',	'NPY4R',	'NPY5R',	'NTSR1',	'NTSR2',	'OPRD1',	'OPRK1',	'OPRL1',	'OPRM1',	'OXTR',	'PRLHR',	'PRLR',	'RAMP1',	'RAMP2',	'RAMP3',	'RXFP1',	'RXFP2',	'SCTR',	'SORT1',	'SSTR1',	'SSTR2',	'SSTR3',	'SSTR4',	'SSTR5',	'TACR1',	'TACR2',	'TACR3',	'TRHR',	'UTS2R',	'VIPR1',	'VIPR2')
gene <- c("RLN1", "GHRL", "CALCA", "PMCH", "PDYN", "UTS2B", "VIP", "GRP", "TAC1", "PENK",
"SCG5", "CBLN4", "TAC4", "NPPC", "GNRH1", "EDN3", "EDN1", "CALCB", "NUCB2",
"NPFF", "NMU", "RLN2", "KNG1", "CHGB", "CBLN2", "CARTPT", "POMC", "NPB", "UCN",
"CORT", "AGT", "SCT", "HCRT", "AGRP", "VGF", "KISS1", "GAST", "UCN2", "GAL",
"TAC3", "CRH", "SST", "ADCYAP1", "PNOC", "CBLN1", "CBLN3", "OXT", "NPW", "GNRH2",
"ADM", "PYY", "NPPB", "TRH", "NPPA", "NMB", "NPY", "CHGA", "CCK", "GHRH", "APLN")

######################加载函数################################
library(dplyr)
fpkm_function <- function(data_frame, reverse = FALSE) {
  if (!reverse) {
    # 这里执行正向操作，例如进行log2转换
    transformed_data <- data_frame %>% 
      mutate(across(-project, ~ 2^(.x) - 0.001))
    return(transformed_data)
  } else {
    # 这里执行反向操作，例如反向转换回原始格式
    reverse_transformed_data <- data_frame %>% 
      mutate(across(-project, ~ log2(.x + 0.001)))
    return(reverse_transformed_data)
  }
}
################################################################
# 提取所需基因的数据并计算均值
##转化为fpkm
data_pro_fpkm <- tcga_mrna_clin %>%
  select(project, all_of(gene))
data_fpkm <- fpkm_function(data_pro_fpkm,reverse = FALSE)
##计算均值
data <- data_fpkm %>%
  group_by(project) %>%
  summarize(across(gene, mean, na.rm = TRUE))
data <- fpkm_function(data,reverse = TRUE)

#write.table(data,file = "纯TCGA肿瘤数据(平均数).csv",sep = ",",row.names = FALSE)
write.table(data,file = "纯TCGA肿瘤数据配体(平均数).csv",sep = ",",row.names = FALSE)
