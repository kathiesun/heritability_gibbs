library(DOQTL)


setwd("/nas02/home/k/y/kys6/pomp_do_intensities")
setwd("/Users/kathiesun/Dropbox (ValdarLab)/outputs/h2_input/")
gm_data <- read.csv("GM_info.csv")


## if founder.probs.Rdata not available
condense.model.probs(path = "DO1_redo", write = "founder.probs.Rdata",
model = "additive",cross = "DO")


load("founder.probs.1.Rdata")
load("founder.probs.2.Rdata")

#####
## for DO1 remove DP.DO1.137.f_8697404046_R11C01 and make
## DP.DO1.137.f_8659977009_R08C01 only mouse 137
remove <- grep("137", rownames(model.probs))[2]
rename <- grep("137", rownames(model.probs))[1]
rownames(model.probs)[rename] <- "DP.DO1.137.f"
model.probs <- model.probs[-remove, , ]

pheno <- read.csv("DO_Pheno9.csv")
covar <- readRDS("DO_covar.rds")
pheno_use <- pheno[which(pheno$DO == 2),]
pheno_use <- pheno_use[,-which(apply(pheno_use, 2, function(x) sum(is.na(x))) > nrow(pheno_use)/3)]
covar_use <- covar[which(covar$DO == 2),]

DO = 2
fem <- ifelse(DO == 1, "f", "F")
mal <- ifelse(DO == 1, "m", "M")

rownames(pheno_use) <- paste0("DP.DO", pheno_use$DO, ".",pheno_use$MouseID, ".",ifelse(pheno_use$Sex == 0, fem, mal))
rownames(covar_use) <- paste0("DP.DO", covar_use$DO, ".", covar_use$MouseID,".",ifelse(covar_use$Sex == "Female", fem, mal))
covar_use$Sex <- ifelse(covar_use$Sex == "Female", 0, 1)
covar_use <- data.matrix(covar_use)
if(length(which(is.na(pheno_use$Sex))) > 0) pheno_use <- pheno_use[-which(is.na(pheno_use$Sex)),]

K = DOQTL::kinship.probs(model.probs, snps = gm_data, bychr = TRUE)

#rownames(model.probs) <- do.call("rbind", strsplit(rownames(model.probs),"[.]"))[,3]

phen_col <- c(6:90)
phen_col <- c(6:98)
qtl = DOQTL::scanone(pheno = pheno_use, pheno.col = phen_col, probs = model.probs, K = K,
                     addcovar = covar_use, snps = gm_data)

saveRDS(qtl, "../h2_outputs/doqtl_scans_do2_again.rds")

quit()



