#!setwd("D:/课本/与同学的合作项目/董子凯/正确数据泛癌/正常与异常组织/正常肿瘤多基因对比图片/差异表达/测试")
#!data<-read.table("GSE98793.txt")
#!group_list<-read.table("gene_list.txt")
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
result <- perform_limma_analysis(data_matrixs, sample_info)
write.table(result, "limma_result.txt")


