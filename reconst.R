library(DOQTL)
library(happy.hbrem)
library(tidyr)
library(magrittr)

##### on killdevil
setwd("~/pomp_do_intensities/alleles")

all1 <- c()
all2 <- c()
for(sim in 1:100){
  all1 <- c(all1, readRDS("DO1_alleles_",sim,".rds"))
  all2 <- c(all2, readRDS("DO2_alleles_",sim,".rds"))
}

saveRDS(all1, "all1.rds")
saveRDS(all2, "all2.rds")

##### on desktop
setwd("/Users/kathiesun")

genot <- read.table("pomp_do_intensities/PompMM08_12222012/UNC-Pomp Mouse 12dec2012_FinalReport_dataOnly.txt", 
                    sep="\t", header=T)
colnames(genot)[3:4] <- c("Allele1", "Allele2")
genot[genot=="-"] = NA  

#tmpvec <- c()
#for(i in 1:nrow(genot)){
#  tmpvec[i] <- ifelse(any(is.na(genot[i,"Allele1"]), is.na(genot[i,"Allele2"])), "N",
#                      ifelse(genot[i,"Allele1"] == genot[i,"Allele2"], as.character(genot[i,"Allele1"]), "H") )
#  print(i)
#}
tmpvec <- readRDS("/pomp_do_intensities/alleles/all2.rds")
genot$allele <- tmpvec  
genot$Sample.ID <- as.character(paste(genot$Sample.ID))
genot_tmp <- genot[,c("SNP.Name","Sample.ID","allele")]

genot_mat <- spread(genot_tmp, SNP.Name, allele, drop=F)  #[unique(genot_mat$SNP.Name),]
colnames(genot_mat) <- genot_mat[1,]
genot_mat <- genot_mat[-1,]

genot_mat <- genot_mat[,grep("DO2", rownames(genot_mat))]
sex <- unlist(lapply(head(strsplit(rownames(genot_mat),"-|_")), function(x) x[[4]]))
names(sex) <- rownames(genot_mat)
gen <- rep(length(sex), "DO10")
names(gen) <- rownames(genot_mat)

list(geno=genot_mat, se)


####
genot_UNL <- read.table("pomp_do_intensities/UNL_083112/UNC-UNL Mega Muga 31aug2012_FinalReport_dataOnly.txt", 
                        sep="\t", header=T)
colnames(genot_UNL)[3:4] <- c("Allele1", "Allele2")

tmpvec <- readRDS("/pomp_do_intensities/alleles/all1.rds")
genot_mat$allele <- tmpvec
genot_tmp <- genot_UNL[,c("SNP.Name","Sample.ID","allele")]
genot_UNL_mat <- spread(genot_tmp, Sample.ID, allele)#[unique(genot_mat$SNP.Name),]

rownames(genot_UNL_mat) <- genot_UNL_mat$SNP.Name
genot_UNL_mat <- genot_UNL_mat[,grep("DO1", colnames(genot_UNL_mat))]
sex_UNL <- toupper(unlist(lapply(strsplit(colnames(genot_UNL_mat),"-|_"), function(x) x[[4]])))

#####

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