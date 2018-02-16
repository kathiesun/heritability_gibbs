library(qtl2)
library(tidyverse)
#### read in data ####
setwd("~/pomp_do_intensities")
pr <- readRDS("rqtl2_do12_probs.rds")
#pr1 <- readRDS("rqtl2_do1_probs.rds")
#pr2 <- readRDS("rqtl2_do2_probs.rds")
pheno <- read.csv("DO_Pheno9.csv")
colnames(pheno)[1] <- "MouseID"
covar <- readRDS("DO_covar.rds")


## from geneseek2qtl2.R
data <- read_cross2("DO12_rqtl2.json", quiet=F)
pr <- calc_genoprob(data, error_prob=0.002)

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
  #filter(MouseID %in% covar[which(covar$DO == 2),]$MouseID) %>%
  select(which(colMeans(is.na(.)) < 0.5)) %>%
  select(-one_of(drop)) %>%
  #select(MouseID, CtClr, XREVS56) %>%
  mutate(CtClr = as.numeric(CtClr)) %>%
  rename(id = MouseID) -> pheno2
#, ends_with('56')

rownames(covar2) <- pheno2$id
rownames(pheno2) <- pheno2$id
#covar1 <- data.matrix(covar1[,-1])
pheno2 <- as.matrix(pheno2[,-1])
pheno2_short <- as.matrix(pheno2[,-c(3:5)])

out2 <- scan1(pr, pheno = pheno2, addcovar = covar2)
outperm <- scan1perm(pr, pheno=pheno2, addcovar=covar2, n_perm=1000)

## plot
color <- c("slateblue", "violetred", "seagreen3", "steelblue1")
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out2) # overall maximum LOD score

for(i in 1:2){
  add=T
  if(i==1) add=F
  plot_scan1(out2, gmap_list, lodcolumn=i, col=color[i], ylim=c(0, ymx*1.02), add=add)
}
legend("topleft", lwd=2, col=paste(color)[1:2], colnames(out2)[1:2], bg="gray90")

find_peaks(out, gmap_list, threshold=4, drop=1.5)

######### calculate kinship mat ######### 
pr1 <- lapply(pr, function(x) x[which(rownames(x) %in% pheno1$MouseID), , ])
pr2 <- lapply(pr, function(x) x[which(rownames(x) %in% pheno2$MouseID), , ])
saveRDS(pr1, "rqtl2_do1_probs.rds")
saveRDS(pr2, "rqtl2_do2_probs.rds")
pr <- readRDS("rqtl2_do12_probs.rds")
apr <- readRDS("rqtl2_do12_alleleprob.rds")
kinship <- calc_kinship(pr)
