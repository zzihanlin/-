rm(list = ls())
library(readxl)
library(data.table)
aa<-read_xlsx('/Users/suwa/Desktop/F240612002/runtime/out/0.1.limma/0.Prepare_Data/胆管癌_protein.xlsx')
bb<-read.csv('/Users/suwa/Desktop/F240612002/runtime/out/0.1.limma/0.Prepare_Data/胆管癌clinical.csv')
# cc<-fread('/Users/suwa/Desktop/F240612002/scRNA-seq联合多组学解析组蛋白/0.Prepare_Data/TCGA-LIHC.survival.tsv.gz')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
library(Seurat)

c(
    list.files("/Users/suwa/Desktop/F240612002/my_custom_scripts/", recursive = T, full.names = T, pattern = "\\.R$")
) %>%
    walk(source)
aa<-as.data.frame(aa)
geneid<-substr(aa$Gene, 1, 15)
   library(clusterProfiler)
   library("org.Hs.eg.db") 
   ### 对基因名进行转换
   ##ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GENETYPE, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIPROT
gene_id <- bitr(aa$Gene, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
aa <- merge(aa, gene_id, by.x = "Gene", by.y = "SYMBOL", all.x = TRUE)

dat<-aa
# ##数据统一
library(dplyr)
library(tidyr)

# Step 1: Identify duplicated genes
dup <- dat$Gene[duplicated(dat$Gene)] %>% unique()

# Step 2: Split duplicated and unique genes
dirt  <- dat %>% filter(Gene %in% dup)
clear <- dat %>% filter(!Gene %in% dup)

# Step 3: Identify expression columns
anno_cols <- c("Gene", "ENSEMBL")
expr_cols <- setdiff(colnames(dirt), anno_cols)

# Step 4: Aggregate duplicates (mean for expr, first for annotation)
dirt_avg <- dirt %>%
  group_by(Gene) %>%
  summarise(
    across(all_of(expr_cols), ~ mean(.x, na.rm = TRUE)),
    ENSEMBL = dplyr::first(.data$ENSEMBL),
    .groups = "drop"
  )

# Step 5: Combine with non-duplicated rows
expr <- bind_rows(clear, dirt_avg)

save(expr,file = '/Users/suwa/Desktop/胆管癌/scRNA-seq联合多组学解析组蛋白/0.Prepare_Data/TCGA_expr.Rdata')
rownames(aa) <- make.unique(as.character(aa$Gene))
aa<-aa[,-1]

clin<-data.frame(sample = bb$'Case.Submitter.ID',group = bb$'Disease.Type',OS = bb$'Days.to.Recurrence',OS.time = bb$'Progression.or.Recurrence')
clin <- data.frame(
  sample = bb$Case.Submitter.ID,
  group = rep("primary", nrow(bb)),  # all samples are "primary"
  OS = bb$Days.to.Recurrence,
  OS.time = bb$Progression.or.Recurrence
)
 clin<-clin[clin$group%in%c('Primary Tumor','Recurrent Tumor','Solid Tissue Normal'),]
 clin$group<-ifelse(clin$group=='Cholangiocarcinoma','Cholangiocarcinoma',NA)
group<-data.frame(sample = colnames(aa))
group<-merge(group,clin,by = 'sample',all = T)
group<-group[-104,]
group$group<-ifelse(group$group=='Cholangiocarcinoma','Tumor','Normal')
write.csv(group,'/Users/suwa/Desktop/胆管癌/data/0.Prepare_Data/group.csv')
group<-read.csv('/Users/suwa/Desktop/胆管癌/data/0.Prepare_Data/group.csv')
group<-group[,c(-1,-4,-5)]
source('/Users/suwa/Desktop/胆管癌/RCodes/limma_deg.R')
expr<-aa[,group$sample]
#' @TODO 差异表达分析
#' @title ### 基于limma进行两组间差异表达分析
#' @param od 结果输出路径
#' @param DEG_exp 表达谱，基因在行样本在列
#' @param DEG_pdata 样本分组文件，第一列为样本，第二列为分组
#' sample   group
#' TCGA-75-6207-01A tumor
#' TCGA-78-7160-01A tumor
#' TCGA-49-6743-01A tumor
#' @param controlLabel 对照组的分类标签
#' @param caseLabel 实验组的分类标签
#' @param DEG_FC 等于log(差异倍数阈值)，默认为1
#' @param DEG_P 差异显著性检验P值
#' @param pvalue 是否使用未校正的P，默认NULL是不使用；如果不是NULL，则使用未校正的P
#' @param saveplot 是否生成结果图片，包括火山图和热图，默认为FALSE不生成
#' @param color_fun 热图中样本分组的颜色
#' color_fun1 <- c("#377EB8", "#E41A1C", "#4DAF4A", "#FF7F00", "#984EA3", "#F781BF", "#A65628", "#FFFF33")
#' @return *list*
#' @examples res <- limma_deg(od = out_dir, DEG_exp = x_exp, DEG_pdata = pdata, controlLabel = "normal", caseLabel = "tumor", DEG_FC = 1, DEG_P = 0.05, color_fun = color_fun1)
#' @author *CY*
#'
expr <- aa[, group$sample]           # subset by sample columns
expr <- as.matrix(expr)              # convert to numeric matrix
mode(expr) <- "numeric"              # ensure all values are numeric

# Mark samples with "N" as Normal, others as Tumor
group$group <- ifelse(grepl("N", group$sample), "Normal", "Tumor")

res <- limma_deg(
  od = "/Users/suwa/Desktop/胆管癌/out/01.limma/",
  DEG_exp = expr,
  DEG_pdata = group,
  controlLabel = "Tumor",
  caseLabel = "Normal",
  DEG_FC = 1.5,
  DEG_P = 0.001,
  pvalue = TRUE,           # ← TRUE to filter using p-value (FALSE to skip)
  saveplot = TRUE,
  color_fun = color_fun4
)

gene<-fread('/Users/suwa/Desktop/胆管癌/out/01.limma/SupplementaryTable_Normal_vs_Tumor_nrDEG.txt')
gene<-as.data.frame(gene)
gene1<-gene[gene$Diff=='up',]$V1
gene2<-gene[gene$Diff=='down',]$V1


##富集分析
#' @TODO 富集分析
#' @title ### 功能富集分析
#' @description 调用了子函数`bubble_plot`、`swr`
#' @param genetype 基因类型，用于输出文件的命名
#' @param od 结果输出路径
#' @param genelist 基因列表
#' @param color_fun 字符串向量，代表颜色
#' \>color_fun
#' "#377EB8" "#E41A1C" "#4DAF4A" "#FF7F00" "#984EA3" "#F781BF" "#A65628"
#' @param pAdjustMethod p值的矫正方法,默认为“BH”.one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @return *list*
#' @examples fun_res <- enrich(genetype = "X", genelist = signature, od = out_dir)
#' @author *CY*
#'
enrich (genetype = 'SYMBOL', genelist = gene2, od = '/Users/suwa/Desktop/胆管癌/out/胆管癌/02.enrich/down/', color_fun = color_fun1, organism = "hsa", pAdjustMethod = "BH", w = 9, h = 6)



library(genefilter)
library(GSVA)

# library(GSVAdata)
library(Biobase)
library(stringr)
library(msigdbr)
h <- msigdbr(species = "Homo sapiens", 
            category = "H")
h <- dplyr::select(h, gs_name, gene_symbol) %>% 
 as.data.frame %>%
split(., .$gs_name) %>%
lapply(., function(x)(x$gene_symbol)) 
gs <- lapply(h, unique)
# devtools::install_version("matrixStats", version = "1.1.0")#回到旧版本1.1.0
# BiocManager::valid()#检查包依赖，确保所有依赖包都已经兼容
##按样本分组
##原发组
# library(matrixStats)
group1<-group[group$group=='Tumor',]$sample
group2<-group[group$group=='Normal',]$sample

expr1<-expr[,colnames(expr)%in%group1]
expr2<-expr[,colnames(expr)%in%group2]

gsva_hallmark1 <- gsva(as.matrix(expr1), gs)
gsva_hallmark2 <- gsva(as.matrix(expr2), gs)
nrow(gsva_hallmark1)
nrow(gsva_hallmark2)


# name<-intersect(rownames(gsva_hallmark1),rownames(gsva_hallmark2))
# gsva_hallmark1<-gsva_hallmark1[name,]
# gsva_hallmark2<-gsva_hallmark2[name,]
gsva_hallmark1<-as.data.frame(gsva_hallmark1)
gsva_hallmark2<-as.data.frame(gsva_hallmark2)
gsva_hallmark<-cbind(gsva_hallmark1,gsva_hallmark2)

write.csv(gsva_hallmark,'/Users/suwa/Desktop/胆管癌/out/03.hallmark/gsva_hallmark.csv')
# pk<-intersect(colnames(gsva_hallmark),clin$sample)
# group<-clin[clin$sample%in%pk,]
# gsva_hallmark<-gsva_hallmark[,pk]
limma_deg (od = '/Users/suwa/Desktop/胆管癌/out/03.hallmark/', DEG_exp = gsva_hallmark, DEG_pdata = group, controlLabel = 'Tumor',
                      caseLabel = 'Normal', DEG_FC = 0, DEG_P = NULL, pvalue = 0.05, saveplot = TRUE, color_fun = color_fun4)


##此处正常该有阳性结果，可惜木有，筛选通路时如果通路不一致可以考虑选择特有通路或者时共有通路，建议使用共有通路

##预后
# ad<-median(as.numeric(gsva_hallmark[rownames(gsva_hallmark)%in%'HALLMARK_PEROXISOME',]))##中值

km<-read.table('/Users/suwa/Desktop/胆管癌/data/0.Prepare_Data/胆管癌gene.txt')
km<-km$V1
exp<-expr[km,bb$Case.Submitter.ID]
cell<-bb
# cell<-cell[,c(-2,-5)]
cell<-as.data.frame(cell)
rownames(cell)<-cell$`Case.Submitter.ID`
cell<-cell[,-1]
cell$Gender<-ifelse(cell$Gender=='male',1,0)
cell[cell$'Primary.Diagnosis'=='Adenocarcinoma, NOS',]$'Primary.Diagnosis' <- 0
cell[cell$'Primary.Diagnosis'=='Lymphoepithelial carcinoma',]$'Primary.Diagnosis' <- 1
cell[cell$'Primary.Diagnosis'=='Mucin-producing Adenocarcinoma',]$'Primary.Diagnosis' <- 2
cell[cell$'Primary.Diagnosis'=='Unknown',]$'Primary.Diagnosis' <- NA
head(cell)
cell<-cell[,c(-3,-4,-8,-9,-13)]

cell[cell$'Tumor.Grade'=='G1',]$'Tumor.Grade' <- 1
cell[cell$'Tumor.Grade'=='G2',]$'Tumor.Grade' <- 2
cell[cell$'Tumor.Grade'=='G3',]$'Tumor.Grade' <- 3
cell[cell$'Tumor.Grade'=='GX',]$'Tumor.Grade' <- NA
cell$'Progression.or.Recurrence'<-ifelse(cell$'Progression.or.Recurrence'=='yes',1,0)

types <- sapply(cell, class)
cell$'Primary.Diagnosis'<-as.numeric(cell$'Primary.Diagnosis')
cell$'Tumor.Grade'<-as.numeric(cell$'Tumor.Grade')
cell$Age.at.Diagnosis<-as.numeric(cell$Age.at.Diagnosis)
cell$Days.to.Recurrence<-as.numeric(cell$Days.to.Recurrence)
cell$Days.to.Last.Follow.Up<-as.numeric(cell$Days.to.Last.Follow.Up)
cell$Days.to.Last.Known.Disease.Status<-as.numeric(cell$Days.to.Last.Known.Disease.Status)

cell<-na.omit(cell)

exp<-as.data.frame(t(exp[,rownames(cell)]))
##计算相关性
data<-cor(exp,cell,method="spearman") #pearson积差相关系数，spearman等级相关系数和kendall秩相关系数
round(data, 2)#保留两位小数#相关性计算
##绘图部分
out_dir<-'/Users/suwa/Desktop/胆管癌/out/04.spss/'
PDR_immune <- cell
# 计算基因与免疫细胞浸润之间的相关性
cor <- psych::corr.test(exp, PDR_immune, method = "spearman", adjust = "none")
r <- cor$r %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    pivot_longer(col = -cell, names_to = "gene", values_to = "cor")
FDR <- cor$p.adj %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    pivot_longer(col = -cell, names_to = "gene", values_to = "FDR")
res <- r %>% inner_join(FDR, by = c("cell", "gene"))
# 作图
res$FDR <- ifelse(res$FDR < 0.05, ifelse(res$FDR < 0.01, "**", "*"), "")
p1 <- ggplot(res, aes(cell, gene)) +
    geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
    scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c") +
    geom_text(aes(label = FDR), col = "black", size = 5) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9)
    ) +
    labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "Correlation"))
plotout(p = p1, od = out_dir, h = 5.5, w = 10, name = "cor_gene_and_clin(PDR)")



##预后
source('/Users/suwa/Desktop/胆管癌/RCodes/make_lasso_cox.r')
#' @TODO 进行LASSO Cox分析
#' @title LASSO Cox分析
#' @param exp  一个data.frame，格式如下
#        TCGA-AA-3534-11A TCGA-AZ-6598-11A TCGA-A6-2678-11A
# ADAM9          4.245974         3.822105         3.942743
# CX3CL1         1.795120         3.249270         2.399230
# HIF1A          4.422943         4.357331         4.470260
# GBP2           3.628283         3.957626         3.757863
# CTSO           3.750734         3.739499         3.627357
# ITGAV          3.388016         3.425854         3.256678
#' @param sur 一个data.frame，格式如下，必须包含  sample 列
#             sample type OS OS.time
# 1 TCGA-3M-AB46-01A STAD  0    1765
# 2 TCGA-3M-AB47-01A STAD  1      NA
# 3 TCGA-B7-5816-01A STAD  0     812
# 4 TCGA-B7-5818-01A STAD  0     356
# 5 TCGA-B7-A5TI-01A STAD  0     595
# 6 TCGA-B7-A5TJ-01A STAD  0     335
#' @param status_col 字符串，生存状态列的列名
#' @param time_col  字符串，生存时间列的列名
#' @param seed  Lasso Cox 分析时设置的种子数，默认为123
#' @param color  作图的颜色，默认为 RColorBrewer::brewer.pal(3, "Dark2")
#' @param out_dir 图片的输出路径
#' @param w 图片的宽度
#' @param h 图片的高度
#' @return 一个data.frame，包含基因及其系数，以及 LASSO Cox 分析的图片结果，图片保存在输出目录下
#   symbol        coef
# 1  ADAM9 -0.12113753
# 2 CX3CL1  0.15463420
# 3 CXCL11 -0.06672390
# 4   SPP1  0.03183605
# 5   REG4 -0.01457143
# 6  PRELP  0.04553819
#' @examples lasso_cox_res <- make_lasso_cox(exp = expression, sur = survival, status_col = "OS", time_col = "OS.time",seed = 123, out_dir = "./", w = 10, h = 3)
#' @author *LLK*
dat<-as.data.frame(t(exp))
clin$OS.time<-ifelse(clin$OS.time=='yes',1,0)
clin<-clin[clin$sample%in%colnames(dat),]
clin$OS<-as.numeric(clin$OS)
clin$OS.time<-as.numeric(clin$OS.time)
rownames(clin)<-NULL
colnames(clin)<-c('sample','group','OS.time','OS')
cox<-make_lasso_cox (exp = dat, sur = clin, status_col = 'OS', time_col = 'OS.time', seed = 123, color = RColorBrewer::brewer.pal(3, "Dark2"), out_dir = '/Pub/Users/guany/Project/木姐的项目/紧急处理/胆管癌/05.cox', w = 10, h = 3)



##高低风险组
source('/Users/suwa/Desktop/胆管癌/RCodes/make_riskscore.r')
#' @TODO 根据LASSO Cox得到的基因及其系数计算riskscore
#' @title 计算riskscore
#' @param exp 基因表达谱，一个data.frame，格式如下
#        TCGA-EB-A41A-01 TCGA-FR-A726-01 TCGA-BF-AAP1-01
# OR4F5          -9.9658         -9.9658         -9.9658
# OR4F29         -5.0116         -3.1714         -9.9658
# OR4F16         -5.0116         -3.1714         -9.9658
# SAMD11          1.4388          0.1903         -2.7274
# NOC2L           4.8748          5.8652          5.5031
# KLHL17          1.8801          2.2391          2.0946
#' @param sur 一个data.frame，格式如下，必须包含  sample 列
#             sample type OS OS.time
# 1 TCGA-3M-AB46-01A STAD  0    1765
# 2 TCGA-3M-AB47-01A STAD  1      NA
# 3 TCGA-B7-5816-01A STAD  0     812
# 4 TCGA-B7-5818-01A STAD  0     356
# 5 TCGA-B7-A5TI-01A STAD  0     595
# 6 TCGA-B7-A5TJ-01A STAD  0     335
#' @param gene_coef  一个data.frame，格式如下,必须包含 symbol 和 coef 列
#' #      symbol         coef
# 1    NFE2L3 -0.046726479
# 2    CTHRC1  0.034844278
# 3     UHRF1 -0.116124798
# 4       F2R  0.004666969
# 5     MMP11  0.034762420
# 6      GPX3  0.010331147
#' @param status_col 字符串，生存状态列的列名
#' @param time_col  字符串，生存时间列的列名
#' @return 一个data.frame，格式如下
#' @examples aa <- make_riskscore(exp = exp, sur = survival, gene_coef = coef,time_col = "OS.time", status_col = "OS")
#' @author *LLK*
#'
score<-make_riskscore (exp = dat, sur = clin, gene_coef = cox, time_col = 'OS.time', status_col = 'OS')##97预后基因表达0，被自动删除
exp<-as.data.frame(exp)
exp<-exp[rownames(score),]
exp$riskscore<-score$riskscore
dat<-as.data.frame(t(exp))
##KM曲线(单基因预后)
source('/Users/suwa/Desktop/胆管癌/RCodes/plot_gene_km.r')
#' @TODO  根据基因表达中位数绘制 KM 曲线
#' @title  根据基因表达中位数绘制 KM 曲线
#' @param exp 一个data.frame,如下所示
#        TCGA-AA-3867-01A TCGA-CA-6719-01A TCGA-NH-A50V-01A TCGA-AA-A01C-01A
# ODC1         6.56209806       5.60079406       5.61674550      6.328781509
# OSR1         0.17272445       0.45802669       0.25419464      0.048964370
# SLC4A7       1.54644070       1.91042811       1.71571595      0.981430367
# DOCK3        0.06911680       0.40356742       0.09082721      0.094624044
# TIMD4        0.04953430       0.30271788       0.19442697      0.001338147
# EPO          0.00274997       0.02575776       0.14117007      0.002749970
#' @param sur 一个data.frame，必须包含  sample 列
#             sample type OS OS.time
# 1 TCGA-3M-AB46-01A STAD  0    1765
# 2 TCGA-3M-AB47-01A STAD  1      NA
# 3 TCGA-B7-5816-01A STAD  0     812
# 4 TCGA-B7-5818-01A STAD  0     356
# 5 TCGA-B7-A5TI-01A STAD  0     595
# 6 TCGA-B7-A5TJ-01A STAD  0     335
#' @param status_col 字符串，生存状态列的列名
#' @param time_col  字符串，生存时间列的列名
#' @param out_dir 字符串向量，结果的输出路径
#' @param best_cutoff  逻辑值，是否使用最佳cut_off分组绘制KM曲线
#' @param minprop  数值型，分组的比例，只有使用cut_off进行分组时，参数才有效
#' @param risk.table  逻辑值，是否绘制样本数量，默认为 FALSE
#' @param h 数值型，图片的高度
#' @param w 数值型，图片的宽度
#' @param color  作图的颜色，默认为 RColorBrewer::brewer.pal(3, "Dark2")
#' @return 所有基因的KM曲线，并保存在输出路径下的KM文件夹下
#' @examples plot_gene_km(exp = exp, sur = TCGA_clinical, status_col = "OS", time_col = "OS.time", out_dir = "./", color = RColorBrewer::brewer.pal(8, "Dark2"))
#' @author *LLK*
#'
clin<-clin[clin$sample%in%rownames(score),]
plot_gene_km (exp = dat, sur = clin, status_col = "OS", time_col = "OS.time", out_dir = '/Users/suwa/Desktop/胆管癌/out/05.cox/', best_cutoff = NULL, minprop = 0.3, risk.table = "FALSE", w = 5, h = 5, color = RColorBrewer::brewer.pal(8, "Dark2")) 


cell<-cell[rownames(score),]
##高低评分组间风险差异
pheno<-as.data.frame(t(cell))
#' @TODO 差异表达分析
#' @title ### 基于limma进行两组间差异表达分析
#' @param od 结果输出路径
#' @param DEG_exp 表达谱，基因在行样本在列
#' @param DEG_pdata 样本分组文件，第一列为样本，第二列为分组
#' sample   group
#' TCGA-75-6207-01A tumor
#' TCGA-78-7160-01A tumor
#' TCGA-49-6743-01A tumor
#' @param controlLabel 对照组的分类标签
#' @param caseLabel 实验组的分类标签
#' @param DEG_FC 等于log(差异倍数阈值)，默认为1
#' @param DEG_P 差异显著性检验P值
#' @param pvalue 是否使用未校正的P，默认NULL是不使用；如果不是NULL，则使用未校正的P
#' @param saveplot 是否生成结果图片，包括火山图和热图，默认为FALSE不生成
#' @param color_fun 热图中样本分组的颜色
#' color_fun1 <- c("#377EB8", "#E41A1C", "#4DAF4A", "#FF7F00", "#984EA3", "#F781BF", "#A65628", "#FFFF33")
#' @return *list*
#' @examples res <- limma_deg(od = out_dir, DEG_exp = x_exp, DEG_pdata = pdata, controlLabel = "normal", caseLabel = "tumor", DEG_FC = 1, DEG_P = 0.05, color_fun = color_fun1)
#' @author *CY*

type<-data.frame(sample = rownames(score),group = score$riskgroup)
limma_deg (od = '/Users/suwa/Desktop/胆管癌/out/06.pheno/', DEG_exp = pheno, DEG_pdata = type, controlLabel = 'Low',
                      caseLabel = 'High', DEG_FC = 0, DEG_P = 0.05, pvalue = NULL, saveplot = TRUE, color_fun = color_fun4)

##箱线图
source('/Users/suwa/Desktop/胆管癌/RCodes/characteristics_plot_by_group.R')
#' @TODO 根据样本分组信息，计算样本特征值的分布差异
#' @title ### 根据样本分组信息，计算样本特征值的分布差异
#' @description 分组间特征值的箱线图、特征值的热图、特征值的分组均值热图
#' @param characteristics_score 矩阵或者数据框，第一列是sample，其它列为特征值
#' \> head(characteristics_score)
#' sample   Bcells    CAFs CD4_Tcells CD8_Tcells Endothelial Macrophages  NKcells
#' TCGA-E… 0.00765 0.0267       0.132     0.0493      0.126      0.00572 1.11e- 8
#' TCGA-H… 0.0147  0.0115       0.261     0.0813      0.0579     0.00695 1.01e- 9
#' TCGA-Y… 0.00648 0.0176       0.130     0.0414      0.0472     0.00267 5.15e-10
#' @param Group 样本分组信息，数据框，行名为样本，只提取第一列的信息作为分组列
#' > head(Group)
#'                 Group
#' TCGA-EJ-5519-01A     A
#' TCGA-HC-7211-01A     A
#' TCGA-Y6-A8TL-01A     A
#' TCGA-EJ-5504-01A     A
#' TCGA-HC-8265-01A     A
#' TCGA-VN-A88O-01A     A
#' @param od 向量，结果输出路径
#' @param color_fun 向量，颜色代码，用于分组的颜色
#' @param type 字符串，用于结果文件的命名
#' @param cluster_rows 逻辑值，默认为FALSE，热图的行是否进行聚类
#' @param heatplot_by_scale 逻辑值，热图展示时是否根据特征进行scale，默认TRUE（均值热图不进行scale）
#' @param saveplot 逻辑值，是否绘图，默认T
#' @return NULL
#' @examples characteristics_plot_by_group(characteristics_score = immune_score_res[[x]], Group = Group, od = od, type = x)
#' @author *CY*
clinical<-cell
clinical$sample<-rownames(clinical)
rownames(clinical)<-NULL
order<- c("sample","Gender" ,                                                
"Primary.Diagnosis"    ,            
"Tumor.Grade"        ,               "Age.at.Diagnosis"   ,              
"Days.to.Recurrence"         ,       "Days.to.Last.Follow.Up"  ,         
"Days.to.Last.Known.Disease.Status" ,"Progression.or.Recurrence"        )
clinical <- clinical[,order]
Group<-data.frame(row.names = type$sample,Group = type$group)
characteristics_plot_by_group (characteristics_score = clinical, Group = Group, od = '/Users/suwa/Desktop/胆管癌/out/06.pheno/', color_fun = color_fun4,
                                          type = 'risk', cluster_rows = FALSE, heatplot_by_scale = TRUE, saveplot = TRUE, heatplot = TRUE) 


# ###高低组临床cox分析
# # source('/Users/suwa/Desktop/胆管癌/RCodes/make_lasso_cox.r')
# source('//Users/suwa/Desktop/胆管癌/scRNA-seq联合多组学解析组蛋白/plot_gene_km.r')
# #' @TODO  根据基因表达中位数绘制 KM 曲线
# #' @title  根据基因表达中位数绘制 KM 曲线
# #' @param exp 一个data.frame,如下所示
# #        TCGA-AA-3867-01A TCGA-CA-6719-01A TCGA-NH-A50V-01A TCGA-AA-A01C-01A
# # ODC1         6.56209806       5.60079406       5.61674550      6.328781509
# # OSR1         0.17272445       0.45802669       0.25419464      0.048964370
# # SLC4A7       1.54644070       1.91042811       1.71571595      0.981430367
# # DOCK3        0.06911680       0.40356742       0.09082721      0.094624044
# # TIMD4        0.04953430       0.30271788       0.19442697      0.001338147
# # EPO          0.00274997       0.02575776       0.14117007      0.002749970
# #' @param sur 一个data.frame，必须包含  sample 列
# #             sample type OS OS.time
# # 1 TCGA-3M-AB46-01A STAD  0    1765
# # 2 TCGA-3M-AB47-01A STAD  1      NA
# # 3 TCGA-B7-5816-01A STAD  0     812
# # 4 TCGA-B7-5818-01A STAD  0     356
# # 5 TCGA-B7-A5TI-01A STAD  0     595
# # 6 TCGA-B7-A5TJ-01A STAD  0     335
# #' @param status_col 字符串，生存状态列的列名
# #' @param time_col  字符串，生存时间列的列名
# #' @param out_dir 字符串向量，结果的输出路径
# #' @param best_cutoff  逻辑值，是否使用最佳cut_off分组绘制KM曲线
# #' @param minprop  数值型，分组的比例，只有使用cut_off进行分组时，参数才有效
# #' @param risk.table  逻辑值，是否绘制样本数量，默认为 FALSE
# #' @param h 数值型，图片的高度
# #' @param w 数值型，图片的宽度
# #' @param color  作图的颜色，默认为 RColorBrewer::brewer.pal(3, "Dark2")
# #' @return 所有基因的KM曲线，并保存在输出路径下的KM文件夹下
# #' @examples plot_gene_km(exp = exp, sur = TCGA_clinical, status_col = "OS", time_col = "OS.time", out_dir = "./", color = RColorBrewer::brewer.pal(8, "Dark2"))
# #' @author *LLK*
# #'
# coxdat<-as.data.frame(t(cell))
# rownames(coxdat)<-c("Gender"  ,                                               
#   "Age"  ,                                                  
#  "Recurrence-free survival (month)"  ,                     
#  "Overall survial (month)" ,                               
#  "Recurrence  (1, yes; 0, no)" ,                           
#  "Survial  (1, dead; 0, alive)" ,                          
# "Liver cirrhosis (1, yes; 0, no)" ,                       
# "Tumor number"      ,                                     
#  "Tumor size (cm)"    ,                                    
# "Lymph node metastasis (1, yes; 0, no)"  ,                
# "Tumor thrombus (1, yes; 0, no)"  ,                       
# "Tumour enapsulation (1, complete; 0, no)"  ,             
# "HBsAg  (1, positive; 0, negative)"  ,                    
# "HBcAb (1, positive; 0, negative)"  ,                     
# "PT, prothrombin time (s)"      ,                         
# "TB, total bilirubin (µmol-L)"   ,                        
# "ALB, albumin, (g-L)"    ,                                
#  "ALT, aminoleucine transferase (U-L)"  ,                  
# "γ-GT, γ-glutamyltransferase (U-L)" ,                     
# "Preoperative  AFP（ng-mL）"  ,                           
# "BCLC stage"    ,                                         
# "TNM stage"       ,                                       
# "Tumor purity by HE staining" ,                           
# "Proteomic subtype"    ,                                  
# "mRNA subtype"       ,                                    
# "IF the same subgroup between Protein and mRNA subtyping",
# "AA signature (1,yes; 0, no )" )

# plot_gene_km (exp = coxdat, sur = clin, status_col = "OS", time_col = "OS.time", out_dir = '/Users/suwa/Desktop/胆管癌/紧急处理/scRNA-seq联合多组学解析组蛋白/06.pheno/', best_cutoff = NULL, minprop = 0.3, risk.table = "FALSE", w = 5, h = 5, color = RColorBrewer::brewer.pal(8, "Dark2"),score = score)


#' @TODO 差异表达分析
#' @title ### 基于limma进行两组间差异表达分析
#' @param od 结果输出路径
#' @param DEG_exp 表达谱，基因在行样本在列
#' @param DEG_pdata 样本分组文件，第一列为样本，第二列为分组
#' sample   group
#' TCGA-75-6207-01A tumor
#' TCGA-78-7160-01A tumor
#' TCGA-49-6743-01A tumor
#' @param controlLabel 对照组的分类标签
#' @param caseLabel 实验组的分类标签
#' @param DEG_FC 等于log(差异倍数阈值)，默认为1
#' @param DEG_P 差异显著性检验P值
#' @param pvalue 是否使用未校正的P，默认NULL是不使用；如果不是NULL，则使用未校正的P
#' @param saveplot 是否生成结果图片，包括火山图和热图，默认为FALSE不生成
#' @param color_fun 热图中样本分组的颜色
#' color_fun1 <- c("#377EB8", "#E41A1C", "#4DAF4A", "#FF7F00", "#984EA3", "#F781BF", "#A65628", "#FFFF33")
#' @return *list*
#' @examples res <- limma_deg(od = out_dir, DEG_exp = x_exp, DEG_pdata = pdata, controlLabel = "normal", caseLabel = "tumor", DEG_FC = 1, DEG_P = 0.05, color_fun = color_fun1)
#' @author *CY*

limma_deg (od = '/Users/suwa/Desktop/胆管癌/out/06.pheno/risk_limma/', DEG_exp = expr, DEG_pdata = type, controlLabel = 'Low',
                      caseLabel = 'High', DEG_FC = 1, DEG_P = 0.05, pvalue = NULL, saveplot = TRUE, color_fun = color_fun4)




##预后差异富集
genelist<-fread('/Users/suwa/Desktop/胆管癌/out/06.pheno/risk_limma/SupplementaryTable_High_vs_Low_nrDEG.txt')
genelist<-as.data.frame(genelist)
genelist1<-genelist[genelist$Diff=='up',]$V1
genelist2<-genelist[genelist$Diff=='down',]$V1


##富集分析
#' @TODO 富集分析
#' @title ### 功能富集分析
#' @description 调用了子函数`bubble_plot`、`swr`
#' @param genetype 基因类型，用于输出文件的命名
#' @param od 结果输出路径
#' @param genelist 基因列表
#' @param color_fun 字符串向量，代表颜色
#' \>color_fun
#' "#377EB8" "#E41A1C" "#4DAF4A" "#FF7F00" "#984EA3" "#F781BF" "#A65628"
#' @param pAdjustMethod p值的矫正方法,默认为“BH”.one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @return *list*
#' @examples fun_res <- enrich(genetype = "X", genelist = signature, od = out_dir)
#' @author *CY*
#'
enrich (genetype = 'SYMBOL', genelist = genelist2, od = '/Users/suwa/Desktop/胆管癌/out/06.pheno/down/', color_fun = color_fun1, organism = "hsa", pAdjustMethod = "BH", w = 9, h = 6)



library(genefilter)
library(GSVA)


# library(GSVAdata)
library(Biobase)
library(stringr)
library(msigdbr)
h <- msigdbr(species = "Homo sapiens", 
            category = "H")
h <- dplyr::select(h, gs_name, gene_symbol) %>% 
 as.data.frame %>%
split(., .$gs_name) %>%
lapply(., function(x)(x$gene_symbol)) 
gs <- lapply(h, unique)
# devtools::install_version("matrixStats", version = "1.1.0")#回到旧版本1.1.0
# BiocManager::valid()#检查包依赖，确保所有依赖包都已经兼容
##按样本分组
##原发组
group1<-type[type$group=='Low',]$sample
group2<-type[type$group=='High',]$sample

expr1<-expr[,colnames(expr)%in%group1]
expr2<-expr[,colnames(expr)%in%group2]

gsva_hallmark1 <- gsva(as.matrix(expr1), gs)
gsva_hallmark2 <- gsva(as.matrix(expr2), gs)
nrow(gsva_hallmark1)
nrow(gsva_hallmark2)


# name<-intersect(rownames(gsva_hallmark1),rownames(gsva_hallmark2))
# gsva_hallmark1<-gsva_hallmark1[name,]
# gsva_hallmark2<-gsva_hallmark2[name,]
gsva_hallmark1<-as.data.frame(gsva_hallmark1)
gsva_hallmark2<-as.data.frame(gsva_hallmark2)
gsva_hallmark<-cbind(gsva_hallmark1,gsva_hallmark2)

write.csv(gsva_hallmark,'/Users/suwa/Desktop/胆管癌/out/06.pheno/hallmark/gsva_hallmark.csv')
# pk<-intersect(colnames(gsva_hallmark),clin$sample)
# group<-clin[clin$sample%in%pk,]
# gsva_hallmark<-gsva_hallmark[,pk]
limma_deg (od = '//Users/suwa/Desktop/胆管癌/out/06.pheno/hallmark/', DEG_exp = gsva_hallmark, DEG_pdata = type, controlLabel = 'Low',
                      caseLabel = 'High', DEG_FC = 0, DEG_P = NULL, pvalue = 0.05, saveplot = TRUE, color_fun = color_fun4)
