<<<<<<< HEAD
#library(happy.hbrem)
library(DOQTL)
library(tidyverse)
library(qtl2)

#### THIS IS ALL FOR DOQTL ####
#### alleles can be used for DOQTL and RQTL though
## code to use probs for ROP in DOQTL in doqtl_scans.R


##### on killdevil
## this gets the alleles out of the original data
## NEVER NEED TO DO THIS AGAIN

setwd("~/pomp_do_intensities")
#setwd("/nas/depts/006/valdar-lab/users/sunk/rqtl_do/")
#source("../../h2_src/doqtl.R")
genot2 <- read.table("/nas/depts/006/valdar-lab/users/sunk/PompMM08_12222012/UNC-Pomp Mouse 12dec2012_FinalReport_dataOnly.txt", 
                    sep="\t", header=T)

all1 <- c()
all2 <- c()
for(sim in 1:100){
  all1 <- c(all1, readRDS("DO1_alleles_",sim,".rds"))
  all2 <- c(all2, readRDS("DO2_alleles_",sim,".rds"))
}

saveRDS(all1, "all1.rds")
saveRDS(all2, "all2.rds")

##### on desktop
setwd(file.path("C:/DB Mount","Dropbox\ (ValdarLab)","outputs"))
genot <- read.table("h2_input/UNC-UNL Mega Muga 31aug2012_FinalReport_dataOnly.txt", sep="\t", header=T)


colnames(genot)[3:4] <- c("Allele1", "Allele2")
genot[genot=="-"] = NA  

#tmpvec <- c()
#for(i in 1:nrow(genot)){
#  tmpvec[i] <- ifelse(any(is.na(genot[i,"Allele1"]), is.na(genot[i,"Allele2"])), "N",
#                      ifelse(genot[i,"Allele1"] == genot[i,"Allele2"], as.character(genot[i,"Allele1"]), "H") )
#  print(i)
#}
#tmpvec <- readRDS("alleles2.rds")

genot$allele <- tmpvec  
genot$Sample.ID <- as.character(paste(genot$Sample.ID))
genot_tmp <- genot[,c("SNP.Name","Sample.ID","allele")]

genot_mat <- spread(genot_tmp, SNP.Name, allele, drop=F)  #[unique(genot_mat$SNP.Name),]
#colnames(genot_mat) <- genot_mat[1,]
#genot_mat <- genot_mat[-1,]
rownames(genot_mat) <- genot_mat[,1]
genot_mat <- genot_mat[,-1]
genot1_intm <- genot1_intm[grep("DO2", rownames(genot1_intm)),]

genot_mat <- genot_mat[grep("DO2", rownames(genot_mat)),]
sex <- unlist(lapply(strsplit(rownames(genot_mat),"-|_"), function(x) x[[4]]))
names(sex) <- rownames(genot_mat)
gen <- rep("DO10", length(sex))
names(gen) <- rownames(genot_mat)

data2 <- list(geno=genot_mat, sex=sex, gen=gen)

## actually ran on longleaf with sh
## sbatch -t 10:00 -n 1 --wrap="R CMD BATCH --no-save runFounderProb.R"

data1 <- readRDS("h2_input/data1.rds")


source("calc.genoprob.R")
calc.genoprob.test(data1, output.dir = "h2_output/do1", 
                      plot = FALSE, sampletype="DO", method="allele")


### intensity
genot <- read.table("/nas/depts/006/valdar-lab/users/sunk/PompMM08_12222012/UNC-Pomp Mouse 12dec2012_FinalReport_dataOnly.txt", 
                    sep="\t", header=T)
remove <- union(which(is.na(genot2$Y)), which(is.na(genot2$X)))
genot2_rem <- genot2[-remove,]
genot2_x <- genot2_rem[,c("SNP.Name","Sample.ID","X")]
genot2_y <- genot2_rem[,c("SNP.Name","Sample.ID","Y")]


genot2_xm <- spread(genot2_x, SNP.Name, X, drop=T)
genot2_ym <- spread(genot2_y, SNP.Name, Y, drop=T)
remove <- union(which(apply(genot2_ym, 2, function(x) sum(is.na(x)) > 0)), 
                which(apply(genot2_xm, 2, function(x) sum(is.na(x)) > 0)))

genot2_xm <- genot2_xm[,-remove]
genot2_ym <- genot2_ym[,-remove]

rownames(genot2_xm) <- genot2_xm[,1]
rownames(genot2_ym) <- genot2_ym[,1]
genot2_xm <- genot2_xm[,-1]
genot2_xm <- genot2_xm[grep("DO2", rownames(genot2_xm)),]
genot2_ym <- genot2_ym[,-1]
genot2_ym <- genot2_ym[grep("DO2", rownames(genot2_ym)),]

sex <- unlist(lapply(strsplit(rownames(genot2_xm),"-|_"), function(x) x[[4]]))
names(sex) <- rownames(genot2_xm)
gen <- rep("DO10", length(sex))
names(gen) <- rownames(genot2_xm)



data2int <- list(x=genot2_xm, y=genot2_ym, sex=sex, gen=gen)
saveRDS(data2int, "data2_int.rds")
data2int <- readRDS("data2_int.rds")

DOQTL:::calc.genoprob(data2int, output.dir = "../do2_int", array = "megamuga",
                      plot = FALSE, sampletype="DO", method="intensity")
####
genot_UNL <- read.table("/nas/depts/006/valdar-lab/users/sunk/UNL_083112/UNC-UNL Mega Muga 31aug2012_FinalReport_dataOnly.txt", 
                        sep="\t", header=T)
colnames(genot_UNL)[3:4] <- c("Allele1", "Allele2")

tmpvec <- readRDS("./alleles/data1.rds")
genot_UNL$allele <- tmpvec
genot_UNL$Sample.ID <- as.character(paste(genot_UNL$Sample.ID))
genot_tmp <- genot_UNL[,c("SNP.Name","Sample.ID","allele")]

genot_UNL_mat <- spread(genot_tmp, SNP.Name, allele)#[unique(genot_mat$SNP.Name),]
rownames(genot_UNL_mat) <- genot_UNL_mat[,1]
genot_UNL_mat <- genot_UNL_mat[,-1]

genot_UNL_mat <- genot_UNL_mat[grep("DO1", rownames(genot_UNL_mat)),]
sex_UNL <- toupper(unlist(lapply(strsplit(rownames(genot_UNL_mat),"-|_"), function(x) x[[4]])))

names(sex_UNL) <- rownames(genot_UNL_mat)
gen <- rep("DO10", length(sex_UNL))
names(gen) <- rownames(genot_UNL_mat)

data1 <- list(geno=genot_UNL_mat, sex=sex_UNL, gen=gen)
    #data1 saved as data1.rds in alleles

DOQTL:::calc.genoprob(data1, output.dir = "do1", 
                      plot = FALSE, sampletype="DO", method="allele")
#calc.genoprob.alleles(data1, output.dir = "../do1")

### intensity

genot1 <- read.table("/nas/depts/006/valdar-lab/users/sunk/UNL_083112/UNC-UNL Mega Muga 31aug2012_FinalReport_dataOnly.txt", 
                        sep="\t", header=T)

genot1_x <- genot1[,c("SNP.Name","Sample.ID","X")]
genot1_y <- genot1[,c("SNP.Name","Sample.ID","Y")]
genot1_xm <- spread(genot1_x, SNP.Name, X, drop=F)
genot1_ym <- spread(genot1_y, SNP.Name, Y, drop=F)
rownames(genot1_xm) <- genot1_xm[,1]
rownames(genot1_ym) <- genot1_ym[,1]
genot1_xm <- genot1_xm[,-1]
genot1_xm <- genot1_xm[grep("DO1", rownames(genot1_xm)),]
genot1_ym <- genot1_ym[,-1]
genot1_ym <- genot1_ym[grep("DO1", rownames(genot1_ym)),]

sex <- unlist(lapply(strsplit(rownames(genot1_xm),"-|_"), function(x) x[[4]]))
names(sex) <- rownames(genot1_xm)
gen <- rep("DO10", length(sex))
names(gen) <- rownames(genot1_xm)

data1int <- list(x=genot1_xm, y=genot1_ym, sex=sex, gen=gen)
saveRDS(data1int, "data1_int.rds")
data1int <- readRDS("data1_int.rds")
DOQTL:::calc.genoprob(data1int, output.dir = "../do1_int", 
                      plot = FALSE, sampletype="DO", method="intensity")

## DOQTL scan with alleles 


## DOQTL scan with intensities from Dan Gatti
qtl = scanone(pheno = pheno, pheno.col = "prop.bm.MN.RET", probs = probs, K = K, 
              addcovar = addcovar, snps = muga_snps)


############# UNUSED? IDK
calc.genoprob(data, chr = "all", output.dir = ".", plot = TRUE, 
              array = c("gigamuga", "megamuga", "muga", "other"), 
              sampletype = c("DO", "CC", "DOF1", "other"), method = c("intensity", "allele"), 
              founders, transprobs, snps)
