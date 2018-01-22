library(DOQTL)

setwd("~/pomp_do_intensities/alleles")

data2int <- readRDS("data2_int.rds")

DOQTL:::calc.genoprob(data2int, output.dir = "../do2_int", 
                      plot = FALSE, sampletype="DO", method="intensity")


data1int <- readRDS("data1_int.rds")

DOQTL:::calc.genoprob(data1int, output.dir = "../do1_int", 
                      plot = FALSE, sampletype="DO", method="intensity")
