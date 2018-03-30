library(qtl2)
library(tidyverse)

######
# this code uses RQTL on previously generated founder probabilities to perform ROP
# using alleles (from DOQTL), just like DOQTL
# code for plotting and other downstream stuff in rqtl_plot.R
#



#### read in data ####
setwd("~/pomp_do_intensities")
pr <- readRDS("rqtl2_do12_probs.rds")
#pr1 <- readRDS("rqtl2_do1_probs.rds")
#pr2 <- readRDS("rqtl2_do2_probs.rds")
pheno <- read.csv("DO_Pheno9.csv")
colnames(pheno)[1] <- "MouseID"
covar <- readRDS("DO_covar.rds")

## create pseudomarker mapping object
gmap <- read.csv("GM/GM/GM_info.csv", stringsAsFactors = T)
gmap$marker <- as.character(gmap$marker)
gmap %>% select(marker, chr, cM) %>%
  filter(chr != "Y" & chr != "M") -> gmap
gmap_list <- list()
for(i in 1:length(unique(gmap$chr))){
  chr <- unique(gmap$chr)[i]
  gmap_list[[i]] <- gmap[which(gmap$chr == chr), 3]
  names(gmap_list[[i]]) <- gmap[which(gmap$chr == chr), "marker"]
}
names(gmap_list) <- unique(gmap$chr)
map <- insert_pseudomarkers(gmap_list, step=1)

######### scans ######### 
rownames(covar) <- covar$MouseID
rownames(pheno) <- pheno$MouseID
covar$Sex <- ifelse(covar$Sex == "Female", 1, 0)
drop <- c("DO", "Wheel", "Wheel","Sex", "Diet")

covar %>% 
  filter(DO == 2) %>% 
  rename(id = MouseID) %>%
  select(Sex) %>%
  as.matrix() -> covar2

pheno %>% 
  filter(DO == 2) %>%
  select(which(colMeans(is.na(.)) < 0.25)) %>%
  select(-one_of(drop)) %>%
  mutate(CtClr = as.numeric(CtClr)) %>%
  rename(id = MouseID) -> pheno2

rownames(covar2) <- pheno2$id
rownames(pheno2) <- pheno2$id
pheno2 <- as.matrix(pheno2[,-1])

out2 <- scan1(pr, pheno = pheno2, addcovar = covar2)
#outperm <- scan1perm(pr, pheno=pheno2, addcovar=covar2, n_perm=1000)

saveRDS(out2, "do2_scan1.rds")
#saveRDS(outperm, "do2_scan1perm.rds")
