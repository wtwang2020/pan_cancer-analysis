# 加载RDATA文件
########################################设置分开代码里的数据路径和数据名称#############################
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据/配体")
#load("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据/TCGA_GTEx_pancancer_mrna_pheno.rdata")

tumor_lable <- c("ACC","BLCA","BRCA","CESC","COAD","DLBC","ESCA","GBM","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM","STAD",
                 "TGCT","THCA","THYM","UCEC","UCS")
library(DESeq2)
library(limma)
library(readr)
library(edgeR)
logFC_t = 2
pvalue_t = 0.05
for (i in 1:length(tumor_lable)) {
  if (i != 4) { 
suffix <- tumor_lable[i]
load(paste(suffix, "配体.RDATA", sep = ""))
gene_expression <- get(paste("gene_expression", suffix, sep = ""))
sample_info <- get(paste("sample_info", suffix, sep = ""))
############################函数加载###############################
perform_limma_analysis <- function(data_matrixs, sample_info) {
  # 数据预处理
  group <- sample_info[, 2]
  group <- factor(group,
                  levels = c("TCGA_tumor","GTEx_normal"))
  suppressMessages(library(limma))
  dge <- DGEList(counts=data_matrixs)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group)
  v <- voom(dge,design, normalize="quantile")   #limma包对RNA-SEQ的话要VOOM标准化
  fit <- lmFit(v, design)
  fit= eBayes(fit)
  DEG3 = topTable(fit, coef=2, n=Inf)     #直接提取 不用as.data.frame
  DEG3 = na.omit(DEG3)
  
  k1 = (DEG3$P.Value < pvalue_t)&(DEG3$logFC < -logFC_t);table(k1)
  k2 = (DEG3$P.Value < pvalue_t)&(DEG3$logFC > logFC_t);table(k2)
  DEG3$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
  table(DEG3$change)
  head(DEG3)

  return(DEG3)
}
#调用函数
#result <- perform_limma_analysis(data_matrixs, sample_info)
#write.table(result, "limma_result.txt")
#################################################################
data_lables <- gene_expression[,1]
data_matrixs <- gene_expression[,-1]
row.names(data_matrixs)<- data_lables
?`limma-package`
data_matrixs <- t(data_matrixs)
print(suffix)
result <- perform_limma_analysis(data_matrixs, sample_info)
write.csv(result, paste(suffix, "配体.csv", sep = ","))
  }}
