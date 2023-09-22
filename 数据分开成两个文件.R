library(tidyverse)
library(dplyr)
library(readr)
#手动转化，加载数据
###########################################更改储存路径和基因集#########################################
#setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析/原始数据存储/count数据")
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析2/原始数据存储/count数据/配体")
setwd("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析2/原始数据存储/count数据/配体")
load("D:/课本/与同学的合作项目/董子凯/无敌稳固泛癌分析2/原始数据存储/count数据/TCGA_GTEx_pancancer_mrna_pheno.rdata")
#genes <- c('ADCYAP1R1',	'ADIPOR1',	'ADIPOR2',	'AGTR1',	'AGTR2',	'APLNR',	'AVPR1A',	'AVPR1B',	'AVPR2',	'BDKRB1',	'BDKRB2',	'CALCR',	'CALCRL',	'CCKAR',	'CCKBR',	'CRHR1',	'CRHR2',	'EDNRA',	'EDNRB',	'GALR1',	'GALR2',	'GALR3',	'GCGR',	'GHR',	'GHRHR',	'GHSR',	'GIPR',	'GLP1R',	'GLP2R',	'GNRHR',	'GRPR',	'HCRTR1',	'HCRTR2',	'KISS1R',	'LEPR',	'MC3R',	'MC4R',	'MCHR1',	'MCHR2',	'MLNR',	'NMBR',	'NMUR1',	'NMUR2',	'NPBWR1',	'NPBWR2',	'NPFFR1',	'NPFFR2',	'NPR1',	'NPR2',	'NPR3',	'NPSR1',	'NPY1R',	'NPY2R',	'NPY4R',	'NPY5R',	'NTSR1',	'NTSR2',	'OPRD1',	'OPRK1',	'OPRL1',	'OPRM1',	'OXTR',	'PRLHR',	'PRLR',	'RAMP1',	'RAMP2',	'RAMP3',	'RXFP1',	'RXFP2',	'SCTR',	'SORT1',	'SSTR1',	'SSTR2',	'SSTR3',	'SSTR4',	'SSTR5',	'TACR1',	'TACR2',	'TACR3',	'TRHR',	'UTS2R',	'VIPR1',	'VIPR2')
genes <- c("RLN1", "GHRL", "CALCA", "PMCH", "PDYN", "UTS2B", "VIP", "GRP", "TAC1", "PENK",
          "SCG5", "CBLN4", "TAC4", "NPPC", "GNRH1", "EDN3", "EDN1", "CALCB", "NUCB2",
          "NPFF", "NMU", "RLN2", "KNG1", "CHGB", "CBLN2", "CARTPT", "POMC", "NPB", "UCN",
          "CORT", "AGT", "SCT", "HCRT", "AGRP", "VGF", "KISS1", "GAST", "UCN2", "GAL",
          "TAC3", "CRH", "SST", "ADCYAP1", "PNOC", "CBLN1", "CBLN3", "OXT", "NPW", "GNRH2",
          "ADM", "PYY", "NPPB", "TRH", "NPPA", "NMB", "NPY", "CHGA", "CCK", "GHRH", "APLN")

data_select <- tcga_gtex_mrna_pheno %>%
  select(sample_id,sample_type,project,primary_site,all_of(genes))
#词频统计################去除重复#######################
keywords_to_remove <- c("CHCL", "NHSC", "MESO", "PCPG", "SARC", "UVM")
# 使用逻辑条件筛选删除行
filtered_data <- data_select[!grepl(paste(keywords_to_remove, collapse = "|"), data_select$project), ]
filtered_data <- filtered_data[!grepl(paste("TCGA_normal", collapse = "|"), filtered_data$sample_type), ]
#################################去除第一列中重复的行##############################
first_column <- filtered_data[, 1]
filtered_data <- filtered_data[!duplicated(first_column), ]
filtered_data_test_frequent <- filtered_data%>%
  filter(sample_type == "GTEx_normal")
duplicates <- duplicated(filtered_data_test_frequent$sample_id)
word_freq <- table(duplicates)
print(word_freq)
##########################数据转化为count##################################
#filtered_data <- filtered_data %>% 
#  mutate(across(-project, ~ 2^(.x))
filtered_data_1 <- filtered_data %>% 
  mutate_at(vars(-sample_id, -sample_type, -project, -primary_site), ~ 2^(.x))

filtered_data <- filtered_data_1 %>%
  mutate_at(vars(-sample_id, -sample_type, -project, -primary_site), ~ floor(.x))

#log_transformed_data <- data_select[, 5:87]
# 应用逆转换公式
#original_count <- apply(log_transformed_data, 2, function(col) 2^col - 1)
# 将结果存储回数据框的相应列
#data_select[, 5:87] <- original_count
#####################删除部分肿瘤以及TCGA正常样本###########################
# 要删除的关键词列表
# 打印筛选后的数据框
##################################################################################
####将第4列中Adrenal Gland替换为ACC,Bladder替换为BLCA，Breast替换为BRCA,Cenvix Uteri替换为CESC，Colon替换为 COAD，Blood替换为DLBC
###选择矩阵
tumor_lable <- c("ACC","BLCA","BRCA","CESC","COAD","DLBC","ESCA","GBM","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SKCM","STAD",
                 "TGCT","THCA","THYM","UCEC","UCS")
nomal_lable <- c("Adrenal Gland","Bladder","Breast","Cenvix Uteri","Colon","Blood","Esophagus","Brain","Kidney","Kidney","Kidney","Bone Marrow","Brain","Liver","Lung","Lung",
                 "Ovary","Pancreas","Prostate","Colon","Skin","Stomach","Testis","Thyroid","Blood","Uterus","Uterus")
for (i in 1:length(tumor_lable)) {
data <- filtered_data %>%
  filter(primary_site %in% c(tumor_lable[i],nomal_lable[i]))
####Esophagus替换为ESCA，Brain替换为GBM，Bone Marrow替换为LAML，#################################
############################
##############################数据分开存储################################
# 提取样本标号、样本来源和样本肿瘤列
####paste("sample_info","_ACC") <- filtered_data[, 1:4]
####sample_info <- filtered_data[, 1:4]
# 创建一个新的数据框用于基因表达信息
####gene_expression <- filtered_data[, -c(2:4)]
##########################################################################
# 假设 sample_info 和 gene_expression 是您的数据框
# 保存数据框为RDATA文件
####save(sample_info, gene_expression, file = "sample_data.RDATA")
suffix <- tumor_lable[i]
new_obj_name_1 <- paste("sample_info", suffix, sep = "")
assign(new_obj_name_1, data[, 1:4])
new_obj_name_2 <- paste("gene_expression", suffix, sep = "")
assign(new_obj_name_2, data[, -c(2:4)])
# 假设 sample_info 和 gene_expression 是您的数据框
# 保存数据框为RDATA文件
#save(list = c(new_obj_name_1, new_obj_name_2), file = paste(suffix, ".RDATA", sep = ""))
save(list = c(new_obj_name_1, new_obj_name_2), file = paste(suffix, "配体.RDATA", sep = ""))
}
