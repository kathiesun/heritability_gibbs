library(happy.hbrem)
library(DOQTL)
library(tidyr)
library(magrittr)
#install.packages("qtl2", repos="https://rqtl.org/qtl2cran")
library(qtl2)


##### on killdevil
setwd("~/pomp_do_intensities/alleles")
#source("../../h2_src/doqtl.R")

all1 <- c()
all2 <- c()
for(sim in 1:100){
  all1 <- c(all1, readRDS("DO1_alleles_",sim,".rds"))
  all2 <- c(all2, readRDS("DO2_alleles_",sim,".rds"))
}

saveRDS(all1, "all1.rds")
saveRDS(all2, "all2.rds")

##### on desktop
#setwd("/Users/kathiesun")

genot <- read.table("../PompMM08_12222012/UNC-Pomp Mouse 12dec2012_FinalReport_dataOnly.txt", 
                    sep="\t", header=T)
colnames(genot)[3:4] <- c("Allele1", "Allele2")
genot[genot=="-"] = NA  

#tmpvec <- c()
#for(i in 1:nrow(genot)){
#  tmpvec[i] <- ifelse(any(is.na(genot[i,"Allele1"]), is.na(genot[i,"Allele2"])), "N",
#                      ifelse(genot[i,"Allele1"] == genot[i,"Allele2"], as.character(genot[i,"Allele1"]), "H") )
#  print(i)
#}
tmpvec <- readRDS("alleles2.rds")
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
data2 <- readRDS("data2.rds")
DOQTL:::calc.genoprob(data2, output.dir = "../do2", 
                      plot = FALSE, sampletype="DO", method="allele")


### intensity
genot2_int <- genot[,c("SNP.Name","Sample.ID","X")]
genot2_intm <- spread(genot2_int, SNP.Name, X, drop=F)
rownames(genot2_intm) <- genot2_intm[,1]
genot2_intm <- genot2_intm[,-1]
genot2_intm <- genot2_intm[grep("DO2", rownames(genot2_intm)),]

sex <- unlist(lapply(strsplit(rownames(genot2_intm),"-|_"), function(x) x[[4]]))
names(sex) <- rownames(genot2_intm)
gen <- rep("DO10", length(sex))
names(gen) <- rownames(genot2_intm)

data2int <- list(geno=genot2_intm, sex=sex, gen=gen)
data2int <- readRDS("data2_int.rds")
DOQTL:::calc.genoprob(data2_int, output.dir = "../do2_int", 
                      plot = FALSE, sampletype="DO", method="intensity")

####
genot_UNL <- read.table("/nas/depts/006/valdar-lab/users/sunk/UNL_083112/UNC-UNL Mega Muga 31aug2012_FinalReport_dataOnly.txt", 
                        sep="\t", header=T)
colnames(genot_UNL)[3:4] <- c("Allele1", "Allele2")

tmpvec <- readRDS("alleles1.rds")
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
calc.genoprob.alleles(data1, output.dir = "../do1")

### intensity
genot1_int <- genot_UNL[,c("SNP.Name","Sample.ID","X")]
genot1_intm <- spread(genot1_int, SNP.Name, X, drop=F)
rownames(genot1_intm) <- genot1_intm[,1]
genot1_intm <- genot1_intm[,-1]
genot1_intm <- genot1_intm[grep("DO1", rownames(genot1_intm)),]

sex <- unlist(lapply(strsplit(rownames(genot1_intm),"-|_"), function(x) x[[4]]))
names(sex) <- toupper(rownames(genot1_intm))
gen <- rep("DO10", length(sex))
names(gen) <- rownames(genot1_intm)

data1int <- list(geno=genot1_intm, sex=sex, gen=gen)
data1int <- readRDS("data1_int.rds")
DOQTL:::calc.genoprob(data1int, output.dir = "../do1_int", 
                      plot = FALSE, sampletype="DO", method="intensity")


####################### OLDER CODE ##########################

calc.genoprob(data, chr = "all", output.dir = ".", plot = TRUE, 
              array = c("gigamuga", "megamuga", "muga", "other"), 
              sampletype = c("DO", "CC", "DOF1", "other"), method = c("intensity", "allele"), 
              founders, transprobs, snps)

DO.output.dir <- "/nas02/home/k/y/kys6/DO_mapping/DO_genotypes"
all.files <- list.files(DO.output.dir)
subject.files <- all.files[grepl(all.files, pattern="genotype.probs.Rdata", fixed=TRUE)]

# from Greg

DO.output.dir <- "~/Documents/Praveen_DOQTL/DOQTL_output/"

all.files <- list.files(DO.output.dir)
subject.files <- all.files[grepl(all.files, pattern="genotype.probs.Rdata", fixed=TRUE)]

subjects <- gsub(subject.files, pattern=".genotype.probs.Rdata", replacement="", fixed=TRUE)

convert.DOQTL.to.HAPPY(DOQTL.recon.output.path=DO.output.dir,
                       map=MM_snps, physical_dist.is.Mb=TRUE, 
                       map.locus_name.colname="SNP_ID", map.chr.colname="Chr", map.physical_dist.colname="Mb_NCBI38", map.genetic_dist.colname="cM",
                       HAPPY.output.path="~/DO_cache/",
                       in.log.scale=FALSE,
                       founder.probs.Rdata.available=TRUE, samples=subjects, simple.alleles=LETTERS[1:8])

library(miqtl)

K <- calc.kinship.from.genomecache.with.DOQTL(your.cache.dir, model="additive")
# Probably good to load up later

DO.data <- read.table(your.data.path, header=TRUE, sep="\t")

ROP.scan <- scan.h2lmm(genomecache=your.cache.dir, data=DO.data,
                       formula=y~1+junk, K=K, model="additive", use.multi.impute=FALSE, 
                       use.fix.par=TRUE, pheno.id=your.id.string)

genome.plotter.whole(scan.list=list(ROP=rop.scan))
# MM_snps.Rdata: Possible map data for DO