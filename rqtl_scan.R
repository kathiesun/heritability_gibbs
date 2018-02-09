setwd("/nas/depts/006/valdar-lab/users/sunk/rqtl_do/")
library(qtl2)
library(tidyverse)


######### scans ######### 
pr1 <- readRDS("rqtl2_do1_probs.rds")
pr2 <- readRDS("rqtl2_do2_probs.rds")
pheno <- readRDS("DO_pheno_only.rds")
covar <- readRDS("DO_covar.rds")
covar$Sex <- ifelse(covar$Sex == "Female", 1, 0)
covar <- data.matrix(covar)
rownames(covar) <- covar[,1]
covar1 <- covar[which(covar[,"DO"] == 1),]
covar2 <- covar[which(covar[,"DO"] == 2),]

pheno %>% filter(MouseID %in% rownames(covar1)) -> pheno1
pheno1 <- select(pheno1, MouseID, starts_with('X'))
pheno1 <- data.matrix(pheno1)
rownames(pheno1) <- pheno1[,1]



covar %>% filter(id %in% covar) -> covar1
covar %>% filter(id %in% rownames(pr2[[1]])) -> covar2
pheno %>% filter(id %in% rownames(pr1[[1]])) -> pheno1
pheno %>% filter(id %in% rownames(pr2[[1]])) -> pheno2

pheno1 <- data.frame(pheno1[match(rownames(pr1[[1]]), pheno1$id), ])
pheno1 <- rbind(colnames(pheno1), pheno1)
out <- scan12(pr1, pheno1, addcovar=covar1, cores=8)

######### calculate kinship mat ######### 
pr1 <- lapply(pr, function(x) x[which(rownames(x) %in% pheno1$MouseID), , ])
pr2 <- lapply(pr, function(x) x[which(rownames(x) %in% pheno2$MouseID), , ])
saveRDS(pr1, "rqtl2_do1_probs.rds")
saveRDS(pr2, "rqtl2_do2_probs.rds")
pr <- readRDS("rqtl2_do12_probs.rds")
apr <- readRDS("rqtl2_do12_alleleprob.rds")
kinship <- calc_kinship(pr)
