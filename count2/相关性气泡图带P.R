##################################data里的基因集和路径，excel大法#################################################
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/图片区/神经肿瘤项目")
# data example
library(dplyr)
library(tidyverse)
library(stringr)
library(ggpubr)
library(tidyr)
library(readr)

#data <- read_csv("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据/差异表达数据.csv", 
#                 col_types = cols(
#                   cancer_type = col_factor(levels = c("ACC","BLCA","BRCA","CESC","COAD","DLBC","ESCA","GBM","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM","STAD",
#                                                       "TGCT","THCA","THYM","UCEC","UCS")),
#                   X = col_factor(levels =  c('ADCYAP1R1',	'ADIPOR1',	'ADIPOR2',	'AGTR1',	'AGTR2',	'APLNR',	'AVPR1A',	'AVPR1B',	'AVPR2',	'BDKRB1',	'BDKRB2',	'CALCR',	'CALCRL',	'CCKAR',	'CCKBR',	'CRHR1',	'CRHR2',	'EDNRA',	'EDNRB',	'GALR1',	'GALR2',	'GALR3',	'GCGR',	'GHR',	'GHRHR',	'GHSR',	'GIPR',	'GLP1R',	'GLP2R',	'GNRHR',	'GRPR',	'HCRTR1',	'HCRTR2',	'KISS1R',	'LEPR',	'MC3R',	'MC4R',	'MCHR1',	'MCHR2',	'MLNR',	'NMBR',	'NMUR1',	'NMUR2',	'NPBWR1',	'NPBWR2',	'NPFFR1',	'NPFFR2',	'NPR1',	'NPR2',	'NPR3',	'NPSR1',	'NPY1R',	'NPY2R',	'NPY4R',	'NPY5R',	'NTSR1',	'NTSR2',	'OPRD1',	'OPRK1',	'OPRL1',	'OPRM1',	'OXTR',	'PRLHR',	'PRLR',	'RAMP1',	'RAMP2',	'RAMP3',	'RXFP1',	'RXFP2',	'SCTR',	'SORT1',	'SSTR1',	'SSTR2',	'SSTR3',	'SSTR4',	'SSTR5',	'TACR1',	'TACR2',	'TACR3',	'TRHR',	'UTS2R',	'VIPR1',	'VIPR2'))
#                 )
#)
data <- read_csv("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据/配体/配体差异表达数据.csv", 
                 col_types = cols(
                   cancer_type = col_factor(levels = c("ACC","BLCA","BRCA","CESC","COAD","DLBC","ESCA","GBM","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM","STAD",
                                                       "TGCT","THCA","THYM","UCEC","UCS")),
                   X = col_factor(levels =  genes <- c("RLN1", "GHRL", "CALCA", "PMCH", "PDYN", "UTS2B", "VIP", "GRP", "TAC1", "PENK",
                                                       "SCG5", "CBLN4", "TAC4", "NPPC", "GNRH1", "EDN3", "EDN1", "CALCB", "NUCB2",
                                                       "NPFF", "NMU", "RLN2", "KNG1", "CHGB", "CBLN2", "CARTPT", "POMC", "NPB", "UCN",
                                                       "CORT", "AGT", "SCT", "HCRT", "AGRP", "VGF", "KISS1", "GAST", "UCN2", "GAL",
                                                       "TAC3", "CRH", "SST", "ADCYAP1", "PNOC", "CBLN1", "CBLN3", "OXT", "NPW", "GNRH2",
                                                       "ADM", "PYY", "NPPB", "TRH", "NPPA", "NMB", "NPY", "CHGA", "CCK", "GHRH", "APLN")
                   ))
)
###########################手动数据预处理################################
############################excel大法###################################
###########################################################################
data <- data[ ,-1]
colnames(data)[1] <- "type"
colnames(data)[2] <- "gene"
colnames(data)[3] <- "Log2FC"
colnames(data)[6] <- "pvalue"
data <- data %>%
  arrange(type, gene)
library(ggplot2)
p1 <- ggplot(data,aes(x=gene ,y= type)) + 
  geom_point(aes(size=-log10(pvalue), fill=Log2FC),
             shape=21,
             color="black") +
  scale_fill_gradient2(name = 'Log2FC\n(Expression)',
                       limit = c(-3.001,3.001),
                       breaks = c(-3,-2,-1,0,1,2,3,4),
                       low='#444283',
                       high='#943934', 
                       mid="white", 
                       midpoint = 0)+
  scale_size_continuous(name = '-Log10 qvalue',
                        limit = c(1.3,7),
                        breaks = c(1,3,5,7))+
  #scale_size_continuous(name = '-Log10 qvalue',
  #                      limit = c(-0.001,3.1),
  #                      breaks = c(0,1,2,3))+
  geom_hline(yintercept=c(5.5, 10.5))+
  labs(x=NULL,
       y=NULL,
       title ="pan-cancer analysis")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 5))+
  scale_x_discrete(labels = function(x) substring(x, 1, 10)) 
  #scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 4)])
p1 
ggsave('bubble_heatmep.pdf',width = 10,height = 5)
