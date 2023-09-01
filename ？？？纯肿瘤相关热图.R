#加载数据
#删除行标签
rownames(data_matrix) <- NULL
#需要的是只有列标签其他都是数值的矩阵
setwd("D:/课本/与同学的合作项目/董子凯/正确数据泛癌")
load("D:/课本/与同学的合作项目/董子凯/正确数据泛癌/相关矩阵.RData")
#install.packages("corrplot")
library(corrplot)
corr <- cor(data_matrix)
#save(data_matrix, file = "相关矩阵.RData")
corrl <- cor.mtest(data_matrix)
corrplot(corr,tl.col = "black", order = "hclust", p.mat = corrl$p,insig = "blank",tl.cex =0.3)
?corrplot
