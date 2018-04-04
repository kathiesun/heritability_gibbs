library(miqtl)
setwd("/nas/depts/006/valdar-lab/users/sunk")
MM_snps <- read.table("DO_mapping/MM_snps.txt")
DO.output.dir <- c("/nas/depts/006/valdar-lab/users/sunk/pomp_do_intensities/DO1_redo",
                   "/nas/depts/006/valdar-lab/users/sunk/pomp_do_intensities/DO2_redo",
                   "/nas02/home/k/y/kys6/DO_mapping/DO_probs_from_gatti")
HAPPY.output.dir <- c("/nas/depts/006/valdar-lab/users/sunk/DO_mapping/miqtl_caches/DO1_cache",
                      "/nas/depts/006/valdar-lab/users/sunk/DO_mapping/miqtl_caches/DO2_cache",
                      "/nas/depts/006/valdar-lab/users/sunk/DO_mapping/miqtl_caches/DO_from_gatti_cache")

for(i in 2:3){
  all.files <- list.files(DO.output.dir[i])
  subject.files <- all.files[grepl(all.files, pattern="genotype.probs.Rdata", fixed=TRUE)]
  subjects <- gsub(subject.files, pattern=".genotype.probs.Rdata", replacement="", fixed=TRUE)
  
  convert.DOQTL.to.HAPPY(DOQTL.recon.output.path=DO.output.dir[i],
                         map=MM_snps, HAPPY.output.path=HAPPY.output.dir[i])
}

convert.DOQTL.to.HAPPY(DOQTL.recon.output.path="/nas/depts/006/valdar-lab/users/sunk/pomp_do_intensities/DO1_redo",
                       map=MM_snps, 
                       HAPPY.output.path="/nas/depts/006/valdar-lab/users/sunk/DO_mapping/miqtl_caches/DO1_cache")






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
