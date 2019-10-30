library(qtl2)

setwd("~/h2_src/")
h2_data <- file.path("C:/DB Mount","Dropbox\ (ValdarLab)","outputs", "h2_outputs")
do2_scan1 <- readRDS(paste0(h2_data, "/do2_scan1.rds"))

## create pseudomarker mapping object (code from rqtl_kildev.R)
gmap <- read.csv(paste0(h2_data, "/../../GM/GM_info.csv"), stringsAsFactors = T)
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
use <- which(colnames(do2_scan1) %in% c("CtClr", "MassPre", "RER_mean", "Feed_mean"))

## plot RQTL
color <- c("slateblue", "violetred", "seagreen3", "steelblue1")
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- max(do2_scan1) # overall maximum LOD score

for(j in 1:ceiling(ncol(do2_scan1)/4)){
  if (j == ceiling(ncol(do2_scan1)/4)) use = ((j-1)*4+1):ncol(do2_scan1)
  else use = c(1:4) + 4*(j-1)
  
  for(i in 1:4){
    add=T
    if(i==1) add=F
    plot_scan1(do2_scan1, gmap_list, lodcolumn=use[i], col=color[i], ylim=c(0, ymx*1.02), add=add)
  }
  
  legend("topleft", lwd=2, col=paste(color)[1:i], colnames(do2_scan1)[use], bg="gray90")
  
}
pks <- find_peaks(do2_scan1, gmap_list, threshold=4, drop=1.5)
hiPks <- pks[which(pks$lod > 30), ]
phenKeep <- unique(hiPks$lodcolumn)[1:4]
use <- which(colnames(do2_scan1) %in% phenKeep)

for(i in 1:length(phenKeep)){
  add=T
  if(i==1) add=F
  plot_scan1(do2_scan1, gmap_list, lodcolumn=use[i], col=color[i], ylim=c(0, ymx*1.02), add=add)
  legend("topleft", lwd=2, col=paste(color)[1:i], colnames(do2_scan1)[use], bg="gray90")
}


######   DOQTL   #######

# do on killdevil

load("do2/founder.probs.Rdata")
pheno <- read.csv("DO_Pheno9.csv")
covar <- readRDS("DO_covar.rds")
pheno2 <- pheno[which(pheno$DO == 2),]
pheno2 <- pheno2[,-which(apply(pheno2, 2, function(x) sum(is.na(x))) > nrow(pheno2)/3)]
covar2 <- covar[which(covar$DO == 2),]

rownames(pheno2) <- paste0("DP.DO", pheno2$DO, ".",pheno2$MouseID, ".",ifelse(pheno2$Sex == 0, "F", "M"))
rownames(covar2) <- paste0("DP.DO", covar2$DO, ".", covar2$MouseID,".",ifelse(covar2$Sex == "Female", "F", "M"))

gmap <- read.csv("GM/GM/GM_info.csv", stringsAsFactors = T)
K = DOQTL::kinship.probs(model.probs, snps = gmap, bychr = TRUE)

#rownames(model.probs) <- do.call("rbind", strsplit(rownames(model.probs),"[.]"))[,3]
qtl = DOQTL::scanone(pheno = pheno2, pheno.col = colnames(pheno2)[6:(dim(pheno2)[2]-1)], probs = model.probs, K = K, 
                     addcovar = covar2, snps = gmap)

save(qtl, "doqtl_out/doqtl_out.rds")



######### MIQTL ###########
# from Greg

DO.output.dir <- "do2"

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




