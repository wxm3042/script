# -------------------------------------
# Date: Fri Oct  7 10:29:36 2022
# Script: 
# Author: WXM
# Purpose: 
# Notes: 
#
# Copyright(c) Corporation Name
# -------------------------------------
library(tidyverse)
library(GENIE3)
library(conflicted)
set.seed(12345)
options(stringsAsFactors = F)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
# detach(unload='conflicted')
# input -------------------------------------------------------------------

tf_path <- 'output/6.interactome/GO_tf/ARACNe/tf.txt'
norm_exp_path <- 'output/2.batchcorrect/norm_gene_count.limma_RBE.txt'

outpath <- 'output/6.interactome/GO_tf/GENIE3/'

tf_list <- read.delim(
  file = tf_path,
  sep = '\t',
  quote = '',
  header = F,
  check.names = F
) %>% nth(1)
norm_exp <- read.delim(
  file = norm_exp_path,
  sep = '\t',
  quote = '',
  header = T,
  check.names = F
) %>% column_to_rownames('gene') %>% 
  as.matrix()

# process -----------------------------------------------------------------

use_tf <- intersect(tf_list, rownames(norm_exp))

weightMat <- GENIE3(norm_exp, regulators=use_tf, nCores=24, verbose=TRUE)

linkList <- getLinkList(weightMat)

write.table(
  linkList,
  file = paste0(outpath, 'network.txt'),
  sep = '\t',
  quote = F,
  row.names = F
)
