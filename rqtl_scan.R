setwd("/nas/depts/006/valdar-lab/users/sunk/rqtl_do/")
library(qtl2)
library(tidyverse)

## read in data
pr <- readRDS("rqtl2_do12_probs.rds")
#pr1 <- readRDS("rqtl2_do1_probs.rds")
#pr2 <- readRDS("rqtl2_do2_probs.rds")
pheno <- readRDS("DO_pheno_only.rds")
covar <- readRDS("DO_covar.rds")

## create pseudomarker mapping object
gmap <- read.csv("GM/GM/GM_info.csv", stringsAsFactors = T)
gmap$marker <- as.character(gmap$marker)
gmap %>% select(marker, chr, cM) -> gmap
gmap_list <- list()
for(i in 1:length(unique(gmap$chr))){
  gmap_list[[i]] <- gmap[which(gmap$chr == 1), 3]
  names(gmap_list[[i]]) <- gmap[which(gmap$chr == 1), "marker"]
}
names(gmap_list) <- unique(gmap$chr)
map <- insert_pseudomarkers(gmap_list, step=1)



######### scans ######### 
rownames(covar) <- covar$MouseID
rownames(pheno) <- pheno$MouseID
covar$Sex <- ifelse(covar$Sex == "Female", 1, 0)

covar %>% 
  filter(DO == 1) %>% 
  select(MouseID, Sex) -> covar1

pheno %>% filter(MouseID %in% covar1$MouseID) %>%
  select(MouseID, ends_with('56')) -> pheno1

rownames(covar1) <- covar1$MouseID
rownames(pheno1) <- pheno1$MouseID
covar1 <- data.matrix(covar1[,-1])
pheno1 <- data.matrix(pheno1[,-1])

out <- scan1(pr, pheno1)

## plot
color <- c("slateblue", "violetred", "seagreen3", "steelblue1")
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out) # overall maximum LOD score

for(i in 1:ncol(out)){
  add=T
  if(i==1) add=F
  plot(out, map, lodcolumn=1, col=color[i], ylim=c(0, ymx*1.02), add=add)
}
legend("topleft", lwd=2, col=paste(color), colnames(out), bg="gray90")

find_peaks(out, map, threshold=4, drop=1.5)

######### calculate kinship mat ######### 
pr1 <- lapply(pr, function(x) x[which(rownames(x) %in% pheno1$MouseID), , ])
pr2 <- lapply(pr, function(x) x[which(rownames(x) %in% pheno2$MouseID), , ])
saveRDS(pr1, "rqtl2_do1_probs.rds")
saveRDS(pr2, "rqtl2_do2_probs.rds")
pr <- readRDS("rqtl2_do12_probs.rds")
apr <- readRDS("rqtl2_do12_alleleprob.rds")
kinship <- calc_kinship(pr)
