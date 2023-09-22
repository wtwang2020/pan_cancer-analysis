
########################################设置读入路径，更改输出名称,更改83,60次数（基因集）#############################
library(readr)
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据/配体")
# 设置文件路径和文件名
tumor_label <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PRAD", "READ", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS")
file_paths <- character(0)  # 初始化一个空的字符向量，用于存储文件名
for (i in 1:length(tumor_label)) {
  if (i != 4) {  # 排除第四个元素
    file_name <- paste(tumor_label[i], ",配体.csv", sep = "")  # 生成文件名
    file_paths <- c(file_paths, file_name)  # 将文件名添加到向量中
  }
}
print(file_paths)
merged_data <- data.frame()
# 创建一个空的列表来存储处理后的数据框
data_list <- list()

for (file_path in file_paths) {
  # 读取当前文件
  current_data <- read.csv(file_path,sep = ",",header = TRUE)  # 根据文件格式进行调整
  
  # 去掉文件名中的 .txt 扩展名
  file_name <- gsub("\\.csv$", "", basename(file_path))
  
  # 将文件名重复83次作为新的列
  #new_column <- rep(file_name, times = 83)
  new_column <- rep(file_name, times = 60)
  
  # 在数据框中插入新列
  current_data <- cbind(cancer_type = new_column, current_data)
  
  # 将处理后的数据框添加到列表中
  data_list[[length(data_list) + 1]] <- current_data
}
# 合并所有处理后的数据框为一个大的数据框
final_data <- do.call(rbind, data_list)
# 打印最终的数据框
print(final_data)
write.csv(final_data, "配体差异表达数据.csv",sep = ",")
#write.csv(final_data, "差异表达数据.csv",sep = ",")
