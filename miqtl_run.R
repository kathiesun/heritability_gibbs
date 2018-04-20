library(miqtl)
setwd("/nas/depts/006/valdar-lab/users/sunk/DO_mapping")
MM_snps <- read.table("MM_snps.txt")
DO.output.dir <- c("../pomp_do_intensities/DO1_redo",
                   "../pomp_do_intensities/DO2_redo",
                   "DO_probs_from_gatti")
HAPPY.output.dir <- c("miqtl_caches/DO1_cache",
                      "miqtl_caches/DO2_cache",
                      "DO_from_gatti_cache")

for(i in 1:3){
  all.files <- list.files(DO.output.dir[i])
  subject.files <- all.files[grepl(all.files, pattern="genotype.probs.Rdata", fixed=TRUE)]
  subjects <- gsub(subject.files, pattern=".genotype.probs.Rdata", replacement="", fixed=TRUE)
  
  convert.DOQTL.to.HAPPY(DOQTL.recon.output.path=DO.output.dir[i],
                         map=MM_snps, HAPPY.output.path=HAPPY.output.dir[i])
}

K <- calc.kinship.from.genomecache.with.DOQTL(your.cache.dir, model="additive")
# Probably good to load up later
caches <- c("miqtl_caches/DO1_cache", "miqtl_caches/DO1_cache", "miqtl_caches/DO_from_gatti_cache")
for(i in 1:3){
  K <- calc.kinship.from.genomecache.with.DOQTL(caches[i], model="additive")
  saveRDS(K, paste0(caches[i],"/kinship.from.DOQTL.rds"))
}

########
library(miqtl)
setwd("/nas/depts/006/valdar-lab/users/sunk/DO_mapping")

DO.data <- read.table("DO_Pheno9.csv", header=TRUE, sep=",")
caches <- c("miqtl_caches/DO1_cache", "miqtl_caches/DO2_cache", "miqtl_caches/DO_from_gatti_cache")
do.dir <- c("DO1_redo", "DO2_redo", "DO_probs_from_gatti")

for(i in 1:3){
  
  K <- readRDS(paste0(caches[i],"/kinship.from.DOQTL.rds")) 
  mouseIDs <- as.numeric(paste(unlist(strsplit(rownames(K),"[.]"))[c(F,F,T,F)]))
  remove <- which(duplicated.default(mouseIDs))
  K <- K[-remove, -remove]
  mouseIDs <- as.numeric(paste(unlist(strsplit(rownames(K),"[.]"))[c(F,F,T,F)]))
  
  use.DO.data <- DO.data[which(DO.data$MouseID %in% mouseIDs),]
  use.DO.data$SUBJECT.NAME <- rownames(K)[match(mouseIDs, use.DO.data$MouseID)]
  keep <- apply(use.DO.data, 2, function(x) ifelse(sum(is.na(x)) < (length(x)/3), T, F))
  use.DO.data <- use.DO.data[,which(keep == T)]

  ROP.scan <- list()

    for(j in 6:ncol(use.DO.data)){
    y <- use.DO.data[,j]
    phen <- colnames(use.DO.data)[j]
    phen.DO.data <- use.DO.data[-which(is.na(y)),]
    K <- K[-which(is.na(y)),-which(is.na(y))]
    y <- y[-which(is.na(y))]
    
    ROP.scan[[phen]] <- scan.h2lmm(genomecache=caches[i], data=phen.DO.data,
                                        formula=y~1+Diet+Sex, K=K, model="additive", 
                                        use.multi.impute=FALSE, use.fix.par=TRUE)
    }
  saveRDS(ROP.scan, paste0(do.dir[i],"/ROPscans.rds"))
}
quit()


#genome.plotter.whole(scan.list=list(ROP=rop.scan))





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

DO.data <- read.table("pomp)do)intensities/DO", header=TRUE, sep="\t")

ROP.scan <- scan.h2lmm(genomecache=your.cache.dir, data=DO.data,
                       formula=y~1+junk, K=K, model="additive", use.multi.impute=FALSE, 
                       use.fix.par=TRUE, pheno.id=your.id.string)

genome.plotter.whole(scan.list=list(ROP=rop.scan))

# MM_snps.Rdata: Possible map data for DO
