#使用easyTCGA获取数据
#https://sfb876.tu-dortmund.de/PublicPublicationFiles/kliewer_lee_2016a.pdf
#清空
rm(list=ls())
gc()
# 安装bioconductor上面的R包
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks")
if(!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment")
if(!require("DESeq2")) BiocManager::install("DESeq2")
if(!require("edgeR")) BiocManager::install("edgeR")
if(!require("limma")) BiocManager::install("limma")
# 安装cran上面的R包
if(!require("survival")) install.packages("survival")
if(!require("broom")) install.packages("broom")
if(!require("devtools")) install.packages("devtools")
if(!require("cli")) install.packages("cli")
#devtools::install_github("ayueme/easyTCGA")
library(easyTCGA)
help(package="easyTCGA")
setwd("F:/data/THCA_TCGA_easyTCGA")
#下载mRNA、lncRNA和临床信息
BRCA<-getmrnaexpr("TCGA-THCA")#原始下载的count, TPM, FPKM 均没有经过log2转化
#下载miRNA
BRCA_miRNA<-getmirnaexpr("TCGA-THCA")
#下载copy number variation data
BRCA_cnv<-getcnv("TCGA-THCA")
#下载masked somatic mutation 体细胞突变
BRCA_snv<-getsnvmaf("TCGA-THCA")
#下载DNA methylation beta value 甲基化数据
getmethybeta("TCGA-THCA")
getclinical("TCGA-THCA")

#从下载目录中打开数据
#差异分析
diff<-diff_analysis(exprset=mrna_expr_counts,#没有经过log2转化
                    project="TCGA-BRCA",
                    save=F)

#批量生存分析
surv<-batch_survival(
  exprset=mrna_expr_counts,
  clin=clin_info,
  is_count = T,
  optimal_cut = TRUE,
  project="TCGA-BRCA",
  save_data = FALSE,
  min_sample_size = 5,
  print_index = TRUE
)

#突变分析：瀑布图
#BiocManager::install("maftools")
library(maftools)
maf<-read.maf(snv,clinicalData=clin_snv)
plotmafSummary(maf)
colnames(clin_snv)
oncoplot(maf=maf,
         clinicalFeatures=c("ajcc_pathologic_stage","vital_status"),
         top=10,
         sortByAnnotation=T
)

#绘制KM曲线
dim(mrna_expr_counts)
set.seed(123)
colnames(clin_info)
clin<-data.frame(time=clin_info$days_to_last_follow_up,
                 event=clin_info$vital_status)
clin$event<-ifelse(clin$event=="Alive",0,1)
plot_KM(exprset=mrna_expr_counts, 
        marker="CHPF", #基因
        clin=clin, 
        optimal_cut = TRUE, 
        return_data = TRUE)

#正常和癌症组织基因表达对比箱线图
rownames(mrna_expr_counts)
plot_gene_paired(exprset=mrna_expr_counts, 
                 marker="CHPF", #基因
                 return_data = TRUE)

#比较组间基因表达差异
set.seed(123)
group=sample(c(0,1),524,replace = T)
plot_gene(exprset=mrna_expr_counts, 
          marker=c("CHPF","MAOA"), 
          group=group, 
          return_data = TRUE)

#单基因差异分析并绘制火山图和热图
#1 加载数据：

load(file = "D:/TCGA/BRCA_easyTCGA/output_mRNA_lncRNA_expr/TCGA-BRCA_mrna_expr_tpm.rdata")
#这个数据是直接从GDC的官网数据中提取出来的，没有经过任何转化，所以我们先进行log2转换。
expr <- log2(mrna_expr_tpm+1)
dim(expr)
#我们这里根据HOPX表达量中位数进行分组，把所有样本分为高表达组和低表达组。
#然后进行差异分析，这里也是用easyTCGA1行代码解决：
sample_group <- ifelse(expr["HOPX",] > median(t(expr["HOPX",])), "high","low")
sample_group <- factor(sample_group, levels = c("low","high"))

library(easyTCGA)
deg_res <- diff_analysis(exprset = expr, 
                         group = sample_group,
                         is_count = F,
                         save = F
)
## => log2 transform not needed
## => Running limma
## => Running wilcoxon test
## => Analysis done.

# 提取limma的结果
deg_limma <- deg_res$deg_limma

head(deg_limma)
#ggplot2绘制火山图
#绘制火山图需要差异分析的结果，我么再增加一列信息展示这个基因是上调、下调还是没意义。
library(dplyr)
library(ggplot2)
library(ggrepel)

deg_limma <- deg_limma %>% 
  mutate(type = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "Up-regulation",
                          logFC < -1 & adj.P.Val < 0.05 ~ "Down-regulation",
                          .default = "None"
  ))
head(deg_limma)
#随便选择几个基因展示一下：
# 要标记的基因
marker <- c("CSTA","AQP3","CD3D","CXCL9","CXCL13","CXCL10","S100A9","S100A8","JCHAIN","KRT5")
#画图
ggplot(deg_limma, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(color=type))+
  scale_color_manual(values = c("blue","gray70","red"),name = NULL)+
  geom_hline(yintercept = -log10(0.05),linetype=2)+
  geom_vline(xintercept = c(-1,1), linetype=2)+
  geom_label_repel(data = subset(deg_limma, genesymbol %in% marker), 
                   aes(label=genesymbol),col="black",alpha = 0.8)+
  theme_bw()+
  theme(legend.position = c(0.15,0.8),
        legend.background = element_blank()
  )

#EnhancedVolcano:一个专门用于绘制火山图的R包，增加了很多好用的功能，比如添加要展示的基因更加方便：
library(EnhancedVolcano)

EnhancedVolcano(deg_limma,
                x = "logFC",
                y = "adj.P.Val",
                lab = rownames(deg_limma),
                selectLab = marker,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                pCutoff = -log10(0.05),
                FCcutoff = 1,
                title = NULL,
                subtitle = NULL
)

#pheatmap绘制热图
#首先是准备热图需要的数据，其实就是表达矩阵的可视化而已。
#我们选择311个差异基因的表达矩阵进行展示。
deg_genes <- deg_limma %>% 
  filter(abs(logFC)>1, adj.P.Val < 0.05) %>% 
  pull(genesymbol)

heat_df <- expr[deg_genes,]

# 自己进行标准化的好处就是可以自定义范围
heat_df <- t(scale(t(heat_df)))
heat_df[heat_df >= 3] <- 3
heat_df[heat_df < -3] <- -3

# 列注释
anno_col <- data.frame(sample_group = sample_group)
rownames(anno_col) <- colnames(heat_df)
#有了数据，画图就很简单：
library(pheatmap)

pheatmap::pheatmap(heat_df,
                   annotation_col = anno_col,
                   show_colnames = F,
                   show_rownames = F
)
#GSEA富集分析可视化
#https://blog.csdn.net/Ayue0616/article/details/132691917?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-1-132691917-blog-132254863.235^v40^pc_relevant_anti_vip&spm=1001.2101.3001.4242.2&utm_relevant_index=4

#TCGA/GTEx泛癌数据1行代码整理
#https://blog.csdn.net/Ayue0616/article/details/131009121?spm=1001.2101.3001.6650.5&utm_medium=distribute.pc_relevant.none-task-blog-2%7Edefault%7EBlogCommendFromBaidu%7ERate-5-131009121-blog-132254863.235%5Ev40%5Epc_relevant_anti_vip&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2%7Edefault%7EBlogCommendFromBaidu%7ERate-5-131009121-blog-132254863.235%5Ev40%5Epc_relevant_anti_vip&utm_relevant_index=8

