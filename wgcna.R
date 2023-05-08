# -------------------------------------
# Date: Thu Jul 14 14:57:01 2022
# Script: 
# Author: WXM
# Purpose: 
# Notes: 
#
# Copyright(c) Corporation Name
# -------------------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)
library(conflicted)
set.seed(12345)
options(stringsAsFactors = F)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
source('src/99.function/KruskalTestFunc.R')
# detach(unload='conflicted')

# function ----------------------------------------------------------------

# datExpr:标准化后表达量数据，行为样本，列为基因
WGCNA_Func <- function(datExpr,output_path,SoftThreshold_cut=0.8,minModuleSize = 100,mergeCutHeight = 0.25){
  library(tidyverse)
  library(WGCNA)
  options(stringsAsFactors = F)
  
  # data input --------------------------------------------------------------
  
  sample_levels <- rownames(datExpr)
  
  # 数据质控 --------------------------------------------------------------------
  
  data_quality <- goodSamplesGenes(datExpr, verbose = 3)
  if(data_quality$allOK){
    print('#######All samples and genes are of good quality!#######')
  }else{
    if(all(data_quality$goodSamples)){
      print('#######All samples are of good quality!#######')
      low_quality_genes_num <- sum(!data_quality$goodGenes)
      print(paste0('#######There are ',low_quality_genes_num,' low quality genes!#######'))
    }else if(all(data_quality$goodGenes)){
      print('#######All genes are of good quality!#######')
      low_quality_samples_num <- sum(!data_quality$goodSamples)
      print(paste0('#######There are ',low_quality_samples_num,' low quality samples!#######'))
    }else{
      low_quality_genes_num <- sum(!data_quality$goodGenes)
      print(paste0('#######There are ',low_quality_genes_num,' low quality genes!#######'))
      low_quality_samples_num <- sum(!data_quality$goodSamples)
      print(paste0('#######There are ',low_quality_samples_num,' low quality samples!#######'))
    }
  }
  datExpr_filter <- datExpr[,data_quality$goodGenes]
  
  # 数据聚类
  sampleTree <-  hclust(dist(datExpr_filter), method = "average")
  pdf(paste0(output_path,'sampleTree.pdf'),width = 6,height = 4)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",cex=0.4)
  dev.off()
  
  # 软阈值 --------------------------------------------------------------------
  
  # 设置线程数
  enableWGCNAThreads(nThreads = 20)
  
  powers <- c(c(1:10),seq(12,30,2))
  sft <-
    pickSoftThreshold(
      datExpr_filter,
      networkType = 'signed',
      powerVector = powers,
      RsquaredCut = SoftThreshold_cut,
      verbose = 5
    )
  nSamples <- nrow(datExpr)
  #画图
  pdf(paste0(output_path,'SoftThreshold.pdf'),width = 10,height = 5)
  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,
       cex=1,
       col="red")
  abline(h=SoftThreshold_cut,col="red")
  #画图
  plot(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       xlab="Soft Threshold (power)",
       ylab="Mean Connectivity", 
       type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], 
       sft$fitIndices[,5], 
       labels=powers, 
       cex=1,
       col="red")
  dev.off()
  
  # 模块鉴定 -------------------------------------------------------------------
  print(sft$powerEstimate)
  net <- blockwiseModules(datExpr_filter,
                          power = sft$powerEstimate,
                          maxBlockSize = 50000,
                          minModuleSize = minModuleSize,
                          TOMType = 'signed',
                          reassignThreshold = 0,
                          mergeCutHeight = mergeCutHeight,
                          numericLabels = T,
                          pamRespectsDendro = F,
                          saveTOMs = T,
                          saveTOMFileBase = paste0(output_path, "blockwiseTOM"),
                          networkType = 'signed',
                          corType = 'pearson',
                          maxPOutliers = 0.1,
                          pearsonFallback = 'individual',
                          verbose = 3)
  saveRDS(net,paste0(output_path,'net.rds'))
  
  # 保存模块基因
  module_gene_mtx <- data.frame(net$colors) %>% 
    rename(color_num=1) %>% 
    rownames_to_column('gene') %>% 
    mutate(module=labels2colors(color_num)) %>% 
    filter(module!='grey') %>% 
    arrange(module)
  write.table(
    module_gene_mtx,
    paste0(output_path,'all_module_gene.txt'),
    sep='\t',
    quote = F,
    row.names = F
  )
  
  # 画图 ----------------------------------------------------------------------
  
  gene_color <- labels2colors(net$colors)
  # 模块聚类图
  pdf(paste0(output_path, 'module.pdf'),
      width = 10,
      height = 5)
  plotDendroAndColors(
    net$dendrograms[[1]],
    gene_color,
    'module colors',
    dendroLabels = F,
    hang = 0.03,
    addGuide = T,
    guideHang = 0.05
  )
  dev.off()
  
  # 模块相关性图
  MEs <- net$MEs
  colnames(MEs) <- colnames(MEs) %>% 
    gsub('^ME','',.) %>% 
    as.numeric() %>% 
    labels2colors(.)
  saveRDS(MEs,paste0(output_path,'MEs.rds'))
  
  pdf(paste0(output_path, 'module_similarity.pdf'),
      width = 7,
      height = 10)
  par(mar = c(6, 12, 3, 1))
  plotEigengeneNetworks(
    MEs,
    "Eigengene adjacency heatmap",
    marDendro = c(0, 5, 5, 5),
    marHeatmap = c(5, 5, 0, 3),
    xLabelsAngle = 90
  )
  dev.off()
  
}

# input -------------------------------------------------------------------

normexp_path <- 'output/2.batchcorrect/norm_gene_count.limma_RBE.txt'
outpath <- 'output/3.WGCNA_allgene/'

normexp <- read.delim(
  file = normexp_path,
  sep = '\t',
  quote = '',
  header = T,
  check.names = F
) %>% 
  select(gene, starts_with('Young'), everything())

group_mtx <- colnames(normexp[,-1]) %>% 
  data.frame(samid=.) %>% 
  mutate(group=gsub('-[0-9]+$','',samid)) %>% 
  column_to_rownames('samid')

kruskal_res <- KruskalTestFunc(normexp, group_mtx)

use_gene <-  kruskal_res %>%
  filter(pval<0.05) %>%
  nth(1)


datexprmtx <- normexp %>% 
  filter(gene %in% use_gene) %>% 
  column_to_rownames('gene') %>% t()

saveRDS(datexprmtx, paste0(outpath, 'datexpr.rds'))

WGCNA_Func(datexprmtx, outpath,mergeCutHeight = 0.4,SoftThreshold_cut=0.8)

# process -----------------------------------------------------------------
# 
# mgene <- read.delim(
#   file = paste0(outpath, 'all_module_gene.txt'),
#   sep = '\t',
#   quote = '',
#   header = T,
#   check.names = F
# )
# 
# plt_list <- list()
# 
# for(m in unique(mgene$module)){
#   subgene <- mgene %>% 
#     filter(module==m) %>% 
#     pull(gene)
#   
#   pdata <- normexp  %>% 
#     filter(gene %in% subgene) %>% 
#     column_to_rownames('gene') %>% 
#     t() %>% scale() %>% t() %>% 
#     data.frame(., check.names = F) %>% 
#     apply(., 2, mean) %>% 
#     data.frame(expr=.) %>% 
#     rownames_to_column('sample_id') %>% 
#     mutate(group=gsub('-[0-9]+$','',sample_id))
#   
#   
#   my_comparisons <- list( c("Young", "Adult"), c("Young", "Aging"))
#   
#   plt_list[[m]] <- ggboxplot(
#     pdata,
#     x='group',
#     y='expr'
#   ) +
#     labs(title=m) +
#     stat_compare_means(comparisons = my_comparisons) +
#     theme(
#       axis.text = element_text(size = 6, color = 'black'),
#       axis.title = element_text(size = 6, color = 'black'),
#       legend.text = element_text(size = 6, color = 'black'),
#       legend.title = element_text(size = 6, color = 'black'),
#       legend.key.size = unit(.1, 'inches')
#     )
# }
# 
# plt_all <- wrap_plots(plt_list,ncol=4)
# 
# ggsave(
#   filename = paste0(outpath, 'module_exp_box.pdf'),
#   plt_all,
#   width = 5,
#   height = 5
# )
