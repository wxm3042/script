# -------------------------------------
# Date: Mon Aug  8 17:15:16 2022
# Script: 
# Author: WXM
# Purpose: 
# Notes: 
#
# Copyright(c) Corporation Name
# -------------------------------------
library(tidyverse)
library(limma)
library(edgeR)
library(sva)
# library(conflicted)
set.seed(12345)
options(stringsAsFactors = F)

# conflict_prefer("filter", "dplyr")
# conflict_prefer("rename", "dplyr")
# detach(unload='conflicted')
# input -------------------------------------------------------------------

raw_path <- 'data/raw_rnaseq_result/gene_symbol_count_matrix.csv'

outpath <- 'output/2.batchcorrect/'

countData <- read.csv(raw_path, check.names=F,row.names = 1)
# process -----------------------------------------------------------------

## 表达logCPM标准化
condition <- sub("-[0-9]+$","",colnames(countData))
d <- DGEList(counts=countData,group=condition) ## 可能受到组别影响
keep <- filterByExpr(d,group=condition,min.count=10,min.total.count=15)
d <- d[keep,,keep.lib.sizes=FALSE]
d <- calcNormFactors(d)
logCPM <- cpm(d,log=TRUE,prior.count=1)

# sva
edata <- logCPM
mod <- model.matrix(~as.factor(condition))
mod0 <- model.matrix(~1,data=as.data.frame(condition))
n.sv <- num.sv(edata,mod,method="leek",seed=12345)
svobj <- sva(edata,mod,mod0,n.sv=n.sv)

# batch correct
design <- model.matrix(~factor(condition))

RBE <- removeBatchEffect(logCPM,
                         covariate = svobj$sv,
                         #batch=as.vector(colnames(logCPM)))
                         design = design)
RBE0 <- cbind(rownames(RBE),RBE)
colnames(RBE0)[1] <- "gene"
write.table(
  RBE0,
  file = paste0(outpath, "norm_gene_count.limma_RBE.txt"),
  sep = "\t",
  quote = F,
  row.names = F
)


