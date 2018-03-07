setwd("/nas02/home/k/y/kys6/pomp_do_intensities")
gm_data <- read.csv("GM/GM/GM_info.csv")
load("do2/founder.probs.Rdata")
pheno <- read.csv("DO_Pheno9.csv")
covar <- readRDS("DO_covar.rds")
pheno2 <- pheno[which(pheno$DO == 2),]
pheno2 <- pheno2[,-which(apply(pheno2, 2, function(x) sum(is.na(x))) > nrow(pheno2)/3)]
covar2 <- covar[which(covar$DO == 2),]

rownames(pheno2) <- paste0("DP.DO", pheno2$DO, ".",pheno2$MouseID, ".",ifelse(pheno2$Sex == 0, "F", "M"))
rownames(covar2) <- paste0("DP.DO", covar2$DO, ".", covar2$MouseID,".",ifelse(covar2$Sex == "Female", "F", "M"))


K = DOQTL::kinship.probs(model.probs, snps = gm_data, bychr = TRUE)

#rownames(model.probs) <- do.call("rbind", strsplit(rownames(model.probs),"[.]"))[,3]
qtl = DOQTL::scanone(pheno = pheno2, pheno.col = "CtClr", probs = model.probs, K = K, 
                     addcovar = covar2, snps = gm_data)