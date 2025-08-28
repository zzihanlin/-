## scRNA数据集 GSE151530_RAW

rm(list=ls())
run_home <- "/Users/suwa/Desktop/cholangiocarcinoma_analysis/"
run_home1 <- '/Users/suwa/Desktop/cholangiocarcinoma_analysis/'
# 导入函数
libSources <- list.files("/Pub/Users/cuiye/RCodes/ProjectCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources, source))
options(bitmapType = "cairo")
pacman::p_load(tidyverse, Seurat, ggplot2, patchwork)

# 设置结果路径
od <- file.path(run_home1, "results/scRNA/GSE151530/1.0.pre_scRNA/")
dir.create(od)
sessionInfo()
# 进入R然后执行下面代码
# devtools::install_version("ggplot2", version = "3.4.4", repos = "http://cran.us.r-project.org")
# devtools::install_version("Matrix", version = "1.5.4")
# devtools::install_github("thomasp85/patchwork")
# BiocManager::install('data.table')
options(stringsAsFactors = F)
options(future.globals.maxSize = 20000 * 1024^2)
library(data.table)
library(COSG)
library(harmony)
library(ggsci)
library(dplyr) 
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(Matrix)
library(Seurat)



path1='/Users/suwa/Desktop/F240612002/runtime/data/GSE151530_RAW/'
list.files(path1) #查看文件明细
sam <- as.data.frame(fread(paste0(path1, "Info.txt.gz")))
head(sam)
head(sam)
#   S_ID Sample               Cell    Type
# 1 S001    H08 AAACCTGAGCTCTCGG-1 T cells
library(readxl)
library(data.table)
aa<-as.data.frame(read_xlsx('/Users/suwa/Desktop/cholangiocarcinoma_analysis/data/scRNA.GSE151530.xlsx',sheet='Sheet1',col_names = F))
dim(aa)
aa=aa[which(aa[,3]=='HCC'),]  # [1] 32  3
head(aa)
#  ...1 ...2 ...3
# 1 GSM4581240 S001  HCC
# 2 GSM4581241 S002  HCC

sam0=sam[which(sam$S_ID %in% aa[,2]),]

Seurat.counts <- Read10X(data.dir = path1) #读取文件
Seurat <- CreateSeuratObject(counts = Seurat.counts, project = "Seurat", min.cells = 3, min.features = 200) #创建Seurat对象，并进行简单的过滤

if (!dir.exists(od)) {
  dir.create(od, recursive = TRUE)
}
saveRDS(Seurat, file = str_glue("{od}/object_GSE151530.rds"))
# > Seurat
# An object of class Seurat 
# 18667 features across 56721 samples within 1 assay 
# Active assay: RNA (18667 features, 0 variable features)
table(Idents(Seurat))
str(Seurat)

labels(Seurat@active.ident)[1:4]
# 1] "AAACCTGAGCTCTCGG-1" "AAACCTGCAAGAAAGG-1" "AAACCTGCAAGTACCT-1"
seurat_obj1=Seurat[, labels(Seurat@active.ident) %in% sam0$'Cell']
str(seurat_obj1)

idents1=data.frame(
  labels(Idents(seurat_obj1)))
colnames(idents1)=c('Cell')

idents2=inner_join(idents1,sam0)
head(idents2)
# Cell S_ID Sample    Type
# AAACCTGAGCTCTCGG-1 S001    H08 T cells

seurat_obj1@meta.data$barcode=idents2$'Cell'
seurat_obj1@meta.data$celltype_ori=idents2$'Type'
seurat_obj1@meta.data$Sample=idents2$'S_ID'
str(seurat_obj1)
saveRDS(seurat_obj1, file = str_glue("{od}/object_GSE151530.hcc.ori.rds"))







# 设置结果路径
od <- file.path(run_home1, "results/scRNA/GSE125449/1.0.pre_scRNA/")
dir.create(od)
path1='/Users/suwa/Desktop/F240612002/runtime/data/GSE125449/'
samples <- list.files(path1)
# 创建一个空的列表来存储Seurat对象
seurat_list <- list()
# 读取每个样本的10x数据并创建Seurat对象
for (sample in samples) {
  # 拼接文件路径
  data.path <- paste0(path1, sample)
  # 读取10x数据，data.dir参数指定存放文件的路径
  seurat_data <- Read10X(data.dir = data.path)
  # 创建Seurat对象，并指定项目名称为样本文件名
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)
  # 将Seurat对象添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

# 打印所有的Seurat对象列表
seurat_list
str(seurat_list[[1]])
# 合并Seurat对象，将所有Seurat对象合并到一个对象中
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)
table(seurat_combined@ meta.data$ orig.ident)
saveRDS(seurat_combined, file = str_glue("{od}/object_GSE125449.rds"))
table(Idents(seurat_combined))
str(seurat_combined)

sam1=as.data.frame(fread('/Users/suwa/Desktop/cholangiocarcinoma_analysis/data/GSE125449/GSE125449_Set1/samples.txt.gz'))
sam2=as.data.frame(fread('/Users/suwa/Desktop/cholangiocarcinoma_analysis/data/GSE125449/GSE125449_Set2/samples.txt.gz'))
head(sam1)
# Sample       Cell Barcode Type
# 1 S02_P01_LCP21 AAACCTGAGGCGTACA-1  CAF
library(readxl)
library(data.table)
aa<-as.data.frame(read_xlsx('/Users/suwa/Desktop/cholangiocarcinoma_analysis/data/scRNA.GSE125449.xlsx',sheet='Sheet1',col_names = F))
dim(aa)
aa=aa[which(aa[,3]=='hcc'),]

sam11=sam1[which(sam1$Sample %in% aa[,2]),]
sam21=sam2[which(sam2$Sample %in% aa[,2]),]

seurat_obj1=seurat_combined[, as.vector(sapply(labels(seurat_combined@active.ident),function(x){substr(x,1+nchar('GSE125449_Set1_'),nchar(x))})) %in% c(sam11$'Cell Barcode',sam21$'Cell Barcode')]
# > seurat_obj1
# An object of class Seurat 
# 21287 features across 3913 samples within 1 assay 
# Active assay: RNA (21287 features, 0 variable features)
str(seurat_obj1)
sam0=rbind(sam11,sam21)
head(sam0)
#   Sample       Cell Barcode Type
# 1 S02_P01_LCP21 AAACCTGAGGCGTACA-1  CAF
# 2 S02_P01_LCP21 AAACGGGAGATCGATA-1  CAF

labels(Idents(seurat_obj1))[1:10]
#  [1] "GSE125449_Set1_AAACCTGAGGCGTACA-1" "GSE125449_Set1_AAACGGGAGATCGATA-1"
#  [3] "GSE125449_Set1_AAAGCAAAGATCGGGT-1" "GSE125449_Set1_AAATGCCGTCTCAACA-1"

idents1=data.frame(
  oriId=labels(Idents(seurat_obj1)),
  newID=as.vector(sapply(labels(Idents(seurat_obj1)),function(x){substr(x,1+nchar('GSE125449_Set1_'),nchar(x))})))
colnames(idents1)=c('oriId','Cell Barcode')

idents2=inner_join(idents1,sam0)
head(idents2)
new.cluster.ids<-idents2$Type

idents2$Type[1:4]
# [1] "CAF" "CAF" "CAF" "CAF"
# https://developer.aliyun.com/article/1256404
seurat_obj1@meta.data$barcode=idents2$'Cell Barcode'
seurat_obj1@meta.data$celltype=idents2$'Type'
seurat_obj1@meta.data$Sample=idents2$'Sample'
str(seurat_obj1)
saveRDS(seurat_obj1, file = str_glue("{od}/object_GSE125449.HCC.celltype.rds"))











