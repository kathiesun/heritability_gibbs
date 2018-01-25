library(DOQTL)
setwd("~/pomp_do_intensities")

genot2 <- read.table("/nas/depts/006/valdar-lab/users/sunk/PompMM08_12222012/UNC-Pomp Mouse 12dec2012_FinalReport_dataOnly.txt", 
                     sep="\t", header=T)
genot1 <- read.table("/nas/depts/006/valdar-lab/users/sunk/UNL_083112/UNC-UNL Mega Muga 31aug2012_FinalReport_dataOnly.txt", 
                     sep="\t", header=T)

## do2
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
#saveRDS(data2int, "data2_int.rds")
#data2int <- readRDS("data2_int.rds")

DOQTL:::calc.genoprob(data2int, output.dir = "alleles/do2_int", array = "megamuga",
                      plot = FALSE, sampletype="DO", method="intensity")

########

## do1

remove <- union(which(is.na(genot1$Y)), which(is.na(genot1$X)))
genot1_rem <- genot1[-remove,]
genot1_x <- genot1_rem[,c("SNP.Name","Sample.ID","X")]
genot1_y <- genot1_rem[,c("SNP.Name","Sample.ID","Y")]

genot1_xm <- spread(genot1_x, SNP.Name, X, drop=T)
genot1_ym <- spread(genot1_y, SNP.Name, Y, drop=T)
remove <- union(which(apply(genot1_ym, 2, function(x) sum(is.na(x)) > 0)), 
                which(apply(genot1_xm, 2, function(x) sum(is.na(x)) > 0)))

genot1_xm <- genot1_xm[,-remove]
genot1_ym <- genot1_ym[,-remove]

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
#saveRDS(data1int, "data1_int.rds")
#data1int <- readRDS("data1_int.rds")

DOQTL:::calc.genoprob(data1int, output.dir = "alleles/do1_int", array = "megamuga",
                      plot = FALSE, sampletype="DO", method="intensity")

