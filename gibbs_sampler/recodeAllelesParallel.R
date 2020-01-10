setwd("~/pomp_do_intensities")
library(cmdline)
sim <- cmdline.integer("sim")

# bsub -J [1-2] R CMD BATCH --vanilla --args --sim=\$LSB_JOBINDEX reconst_par.R
# submit sim 1:100

genot <- read.table("PompMM08_12222012/UNC-Pomp Mouse 12dec2012_FinalReport_dataOnly.txt", 
                    sep="\t", header=T)
colnames(genot)[3:4] <- c("Allele1", "Allele2")
genot[genot=="-"] = NA  

output_dir <- "~/pomp_do_intensities/alleles"

runAll <- function(df,
                   output_dir,
                   perJob=300000,
                   sim) { 
  a = (perJob * (sim-1)) + 1
  b = a+perJob - 1
  if(sim==(floor(nrow(df)/perJob)+1)) b = nrow(df)
  
  mat <- df[a:b,]
  tmpvec <- c()
  for(i in 1:nrow(mat)){
    tmpvec[i] <- ifelse(any(is.na(mat[i,"Allele1"]), is.na(mat[i,"Allele2"])), "N",
                        ifelse(mat[i,"Allele1"] == mat[i,"Allele2"], as.character(mat[i,"Allele1"]), "H") )
  }
  saveRDS(test, paste(output_dir, paste0("DO2_alleles_",sim,".rds"), sep="/"))
  return(tmpvec)
}


test <- runAll(df=genot, 
               output_dir=output_dir,
               sim=sim)


