# -------------------------------------
# Date: Wed Nov 30 16:21:44 2022
# Script: 
# Author: WXM
# Purpose: 
# Notes: 
#
# Copyright(c) Corporation Name
# -------------------------------------
library(tidyverse)
library(ggVennDiagram)
library(ggrepel)

library(conflicted)
set.seed(12345)
options(stringsAsFactors = F)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
# detach(unload='conflicted')
# input -------------------------------------------------------------------

hsa2cfa_path <- '/media/sdd/wxm/database/NCBI/Gene_Transfer/20220713/cfa/hsa2cfa_symbol'
tf_path <- 'output/6.interactome_2/motifs-v10nr_clust-nr.hgnc_tf'
hub_list <- list(
  blue='output/3.WGCNA_allgene/hubgene/bluehub_gene.txt',
  turquoise='output/3.WGCNA_allgene/hubgene/turquoisehub_gene.txt'
)
degree_data_list <- list(
  blue='output/3.WGCNA_allgene/bluecentrality.txt',
  turquoise='output/3.WGCNA_allgene/turquoisecentrality.txt'
)
color_list <- list(
  blue=c('#18499e','#caf0f8'),
  turquoise=c('#61c0ba','#c4fff9')
)

outpath <- 'output/3.WGCNA_allgene/hubgene/'

hsa2cfa <- read.delim(
  file = hsa2cfa_path,
  sep = '\t',
  quote = '',
  header = F,
  check.names = F
) %>% select(3,7) %>% 
  rename(hsa=1, cfa=2)

tf_list <- read.delim(
  file = tf_path,
  sep = '\t',
  quote = '',
  header = F,
  check.names = F
) %>% nth(1)
# process -----------------------------------------------------------------

pdata_list <- list()
for(g in names(hub_list)){
  degree_df <- read.delim(
    file = degree_data_list[[g]],
    sep = '\t',
    quote = '',
    header = T,
    check.names = F
  ) %>% merge(., hsa2cfa, by.x='gene',by.y='cfa') %>% 
    select(hsa, everything(), -gene) %>% 
    rename(gene=1)
  mhubgene <- read.delim(
    file = hub_list[[g]],
    sep = '\t',
    quote = '',
    header = T,
    check.names = F
  ) %>% merge(., hsa2cfa, by.x='gene', by.y='cfa') %>% 
    pull(hsa)
  pdata_list[[g]] <- mhubgene
  hub_tf <- intersect(mhubgene, tf_list)
  hub_tf_degree <- degree_df %>% 
    filter(gene %in% hub_tf)
  write.table(
    hub_tf_degree,
    file = paste0(outpath, g, '_hub_TF.txt'),
    sep = '\t',
    quote = F,
    row.names = F
  )
  pdata <- hub_tf_degree %>% 
    select(-show) %>% 
    arrange(mrank) %>% 
    mutate(rank=1:nrow(.)) %>% 
    mutate(show=ifelse(rank<16, gene, ''))
  plt <- ggplot(pdata, aes(degree, BC, color=rank)) +
    geom_point() +
    geom_label_repel(aes(label=show), size=2, max.overlaps = 100) +
    scale_x_sqrt() +
    scale_color_gradient(low=color_list[[g]][1], high = color_list[[g]][2],) +
    labs(x='Degree', y='BC') +
    theme_bw() +
    theme(
      axis.text = element_text(size = 6, color = 'black'),
      axis.title = element_text(size = 6, color = 'black'),
      legend.text = element_text(size = 6, color = 'black'),
      legend.title = element_text(size = 6, color = 'black'),
      legend.position = 'none',
      panel.grid = element_blank()
    )
  ggsave(
    filename = paste0(outpath, g, '_hub_tf.pdf'),
    plt,
    width = 3,
    height = 3
  )
}
pdata_list[['tf']] <- tf_list
plt <- ggVennDiagram(pdata_list)
ggsave(
  filename = paste0(outpath, 'venn.pdf'),
  plt,
  width = 5,
  height = 5
)
