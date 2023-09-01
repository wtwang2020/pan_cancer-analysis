# 加载RDATA文件
########################################设置分开代码里的数据路径和数据名称#############################
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据/配体")
#load("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据/TCGA_GTEx_pancancer_mrna_pheno.rdata")

tumor_lable <- c("ACC","BLCA","BRCA","CESC","COAD","DLBC","ESCA","GBM","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM","STAD",
                 "TGCT","THCA","THYM","UCEC","UCS")
library(limma)
library(readr)
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
  suppressMessages(library(limma))
  design <- model.matrix(~0 + factor(group))
  colnames(design) <- levels(factor(group))
  rownames(design) <- colnames(data_matrixs)
  contrast.matrix <- makeContrasts("TCGA_tumor-GTEx_normal", levels = design)
  
  # Step 1
  fit <- lmFit(data_matrixs, design)
  
  # Step 2
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Step 3
  tempOutput <- topTable(fit2, coef = 1, n = Inf)
  nrDEG <- na.omit(tempOutput)
  
  return(nrDEG)
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
