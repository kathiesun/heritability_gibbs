library(DOQTL)

##### on killdevil
setwd("~/pomp_do_intensities/alleles")

data2 <- readRDS("data2.rds")
DOQTL:::calc.genoprob(data2, output.dir = "../do2", 
                      plot = FALSE, sampletype="DO", method="allele")

data1 <- readRDS("data1.rds")
DOQTL:::calc.genoprob(data2, output.dir = "../do1", 
                      plot = FALSE, sampletype="DO", method="allele")

quit()
