# -------------------------------------
# Date: Mon Sep 05 14:09:05 2022
# Script: 
# Author: WXM
# Purpose: 
# Notes: 
#
# Copyright(c) Corporation Name
# -------------------------------------
library(tidyverse)
library(fgsea)
library(conflicted)
set.seed(12345)
options(stringsAsFactors = F)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
source('src/99.function/PathwayClusterFunc.R')
source('src/99.function/ReadUnevenFile.R')
source('src/99.function/replotfgseaRes.R')
# detach(unload='conflicted')
# input -------------------------------------------------------------------

deg_path_list <- list(
  # Young_vs_Adult='output/2.batchcorrect/Young_vs_Adult.total_genes_hsa.xls',
  # Adult_vs_Aging='output/2.batchcorrect/Adult_vs_Aging.total_genes_hsa.xls',
  Young_vs_Aging='output/2.batchcorrect/Young_vs_Aging.total_genes_hsa.xls'
)
gmt_path <- '../../database/MSigDB/v7.5.1/c2.cp.reactome.v7.5.1.symbols.gmt'

outpath <- 'output/2.batchcorrect/enrich/fgsea/Reactome/'


gmt <- ReadUnevenFile(gmt_path)

gmt_list <- gmt %>% 
  rename(name=1) %>% 
  select(-2) %>% 
  column_to_rownames('name') %>% 
  apply(., 1, function(x){
    as.character(x) %>% unique(.) %>% setdiff(., '')
  })


# process -----------------------------------------------------------------

for(g in names(deg_path_list)){
  ranks <- read.delim(
    file = deg_path_list[[g]],
    sep = '\t',
    quote = '',
    header = T,
    check.names = F
  ) %>% select(gene, logFC) %>% 
    arrange(logFC) %>% 
    mutate(logFC=-1*logFC) %>% 
    deframe(.)
  
  fgseaRes <- fgsea(
    gmt_list,
    ranks,
    minSize = 5,
    maxSize = 5000
    # nperm = 1000
  )
  
  allgsea <- fgseaRes %>% 
    arrange(NES) %>% 
    mutate(leadinggene=lapply(leadingEdge, function(x){paste0(x, collapse = ';')}) %>% unlist()) %>% 
    mutate(hitnum=lapply(leadingEdge, length) %>% unlist()) %>% 
    filter(hitnum>3) %>% 
    select(-leadingEdge)
  
  write.table(
    allgsea,
    file = paste0(outpath, g, '_allgsea_reactome.xls'),
    sep = '\t',
    quote = F,
    row.names = F
  )
  
  
  siggsea <- allgsea %>% 
    filter(pval<0.05)
  
  write.table(
    siggsea,
    file = paste0(outpath, g, '_siggsea_reactome.xls'),
    sep = '\t',
    quote = F,
    row.names = F
  )
  
}


# cluster

for(g in names(deg_path_list)){
  enrich_res <- read.delim(
    file = paste0(outpath, g, '_siggsea_reactome.xls'),
    sep = '\t',
    quote = '',
    header = T,
    check.names = F
  ) %>% select(pathway, pval, NES, leadinggene) %>% 
    rename(pvalue=2)
  
  pathway_cluster <- PathwayClusterFunc(gmt, enrich_res)
  
  write.table(
    pathway_cluster$net,
    file = paste0(outpath, g, '_net.txt'),
    sep = '\t',
    quote = F,
    row.names = F
  )
  write.table(
    pathway_cluster$node,
    file = paste0(outpath, g, '_node.txt'),
    sep = '\t',
    quote = F,
    row.names = F
  )
  
}




