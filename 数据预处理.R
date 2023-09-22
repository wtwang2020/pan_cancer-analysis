# 加载RDATA文件
########################################设置分开代码里的数据路径和数据名称#############################
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析2/原始数据存储/count数据/配体")
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
#####test
#gene_expression <- gene_expressionACC
#sample_info <- sample_infoACC
#data_lables <- gene_expression[,1]
#data_matrixs <- gene_expression[,-1]
#row.names(data_matrixs)<- data_lables
#?`limma-package`
#data_matrixs <- t(data_matrixs)

############################函数加载###############################
perform_deseq2_analysis <- function(data_matrixs, sample_info) {
  # 数据预处理
  group <- sample_info[,c(1,2)]
  suppressMessages(library(DESeq2))
  #BiocManager::install("apeglm")
  #install.packages("apeglm")
  #suppressMessages(library(apeglm))
  dds <- DESeqDataSetFromMatrix(data_matrixs, group, ~sample_type)
  keep <- rowSums(counts(dds)) >=10
  dds <- dds[keep,]
  ddsDE <- DESeq(dds)
  normCounts <- counts(ddsDE,normalized = T) 
  res <- results(ddsDE, alpha = 0.05)
  summary(res)
  resOrdered <- res[order(res$padj),]
  write.csv(resOrdered,"临时文件.csv")
  decseq_data <- read_csv("临时文件.csv",col_names = TRUE)
  decseq_data[,1]
  decseq_data <- as.data.frame(decseq_data)
  row.names(decseq_data) <- decseq_data[,1]
  return(decseq_data)
}
#调用函数
#result <- perform_limma_analysis(data_matrixs, sample_info)
#write.table(result, "limma_result.txt")
#################################################################
data_lables <- gene_expression[,1]
data_matrixs <- gene_expression[,-1]
row.names(data_matrixs)<- data_lables
data_matrixs <- t(data_matrixs)
print(suffix)
result <- perform_deseq2_analysis(data_matrixs, sample_info)
write.csv(result, paste(suffix, "配体.csv", sep = ","))
  }}
