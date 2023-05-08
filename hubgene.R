# -------------------------------------
# Date: Mon Nov 28 21:36:38 2022
# Script: 
# Author: WXM
# Purpose: 
# Notes: 
#
# Copyright(c) Corporation Name
# -------------------------------------
library(tidyverse)

library(conflicted)
set.seed(12345)
options(stringsAsFactors = F)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
# detach(unload='conflicted')
# input -------------------------------------------------------------------

net_path <- 'output/3.WGCNA_allgene/net.rds'
me_path <- 'output/3.WGCNA_allgene/MEs.rds'
module_gene_path <- 'output/3.WGCNA_allgene/all_module_gene.txt'
dat_path <- 'output/3.WGCNA_allgene/datexpr.rds'

outpath <- 'output/3.WGCNA_allgene/hubgene/'

datexpr <- readRDS(dat_path)
module_gene <- read.delim(
  file = module_gene_path,
  sep = '\t',
  quote = '',
  header = T,
  check.names = F
)
net <- readRDS(net_path)
MEs <- readRDS(me_path)

hsa2cfa <- read.delim(
  file = hsa2cfa_path,
  sep = '\t',
  quote = '',
  header = F,
  check.names = F
) %>% select(3, 7) %>%
  rename(hsa = 1, cfa = 2)

# process -----------------------------------------------------------------

for(m in unique(module_gene$module)){
  subm_genes <- module_gene %>% 
    filter(module == m) %>% 
    pull(gene)
  sub_exp <- datexpr[, subm_genes]
  mmcor <- cor(sub_exp, MEs[,m]) %>% 
    data.frame(., check.names = F) %>% 
    rename(corr=1) %>% 
    rownames_to_column('gene') %>% 
    mutate(core=ifelse(corr>0.7, 'y', 'n'))
  write.table(
    mmcor,
    file = paste0(outpath, m, '_MM.xls'),
    sep = '\t',
    quote = F,
    row.names = F
  )
  hub_gene <- mmcor %>% 
    filter(corr>0.7) %>% 
    pull(gene)
  # write.table(
  #   hub_gene,
  #   file = paste0(outpath, m, 'hub_gene.txt'),
  #   sep = '\t',
  #   quote = F,
  #   row.names = F,
  #   col.names = 'gene'
  # )
}


# plot --------------------------------------------------------------------

hub_gene_list <-list()
for(m in unique(module_gene$module)){
  mm_corr <- read.delim(
    file = paste0(outpath, m, '_MM.xls'),
    sep = '\t',
    quote = '',
    header = T,
    check.names = F
  )
  degree_score <- read.delim(
    file = paste0('output/3.WGCNA_allgene/', m, 'centrality.txt'),
    sep = '\t',
    quote = '',
    header = T,
    check.names = F
  ) %>% mutate(core_d=ifelse(rank<nrow(.)/2, 'y', 'n')) %>% 
    select(-show)
  
  pdata <- merge(mm_corr, degree_score, by='gene') %>% 
    mutate(show=ifelse(core=='y' & core_d=='y', 'y', 'n'))
  write.table(
    pdata,
    file = paste0(outpath, m, '_pdata.xls'),
    sep = '\t',
    quote = F,
    row.names = F
  )
  hub_gene_list[[m]] <- pdata %>% 
    filter(show=='y') %>% 
    mutate(module=m) %>% 
    select(gene, module)
  plt <- ggplot(pdata, aes(corr, degree, color=show)) +
    geom_point() +
    scale_color_manual(values = c(y=m, n='grey')) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 6, color = 'black'),
      axis.title = element_text(size = 6, color = 'black'),
      legend.text = element_text(size = 6, color = 'black'),
      legend.title = element_text(size = 6, color = 'black'),
      legend.key.size = unit(.1, 'inches')
    )
  ggsave(
    filename = paste0(outpath, m, '_hubgene_dot.pdf'),
    plt,
    width = 5,
    height = 5
  )
  
}

all_hub <- bind_rows(hub_gene_list)
write.table(
  all_hub,
  file = paste0(outpath, 'all_hub.txt'),
  sep = '\t',
  quote = F,
  row.names = F
)

for(g in names(hub_gene_list)){
  mhub <- hub_gene_list[[g]] %>% 
    pull(gene)
  write.table(
    mhub,
    file = paste0(outpath, g, 'hub_gene.txt'),
    sep = '\t',
    quote = F,
    row.names = F,
    col.names = 'gene'
  )
  
}

