#加载R包；
library(readr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(superheat)
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/图片区/神经肿瘤项目")
################读取dataframe数据########################
####################更改读取以及输出##########################
#data <- read_csv("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/TCGA纯癌症数据矩阵（基因已经选择）/纯TCGA肿瘤数据(平均数).csv")
data <- read_csv("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/TCGA纯癌症数据矩阵（基因已经选择）/纯TCGA肿瘤数据配体(中位数).csv")
#################################转化为矩阵##############
data_lable <- data[ ,1]
data_lable <- data_lable$project
data <- data[ ,-1]
rownames(data)<- data_lable
data_matrix <- as.matrix(data)
#热图1
compare_heatmap(data_matrix,scale = "row")
#热图2
#pdf('纯肿瘤热图TCGA（平均数）.pdf',height = 5.5,width = 12)
pdf('纯肿瘤热图TCGA配体（中位数）.pdf',height = 5.5,width = 12)
pheatmap(data_matrix,
         color = colorRampPalette(c('#003399', "white", '#CC0033'))(100),
         fontsize_col = 5.5,
         fontsize_row = 8)
dev.off()
?pheatmap
