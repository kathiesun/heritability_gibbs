library(qtl2)

do2_scan <- readRDS("E:/Dropbox\ (ValdarLab)/outputs/h2_outputsdo2_scan1.rds")
use <- which(colnames(do2_scan) %in% c("CtClr", "MassPre", "MassPost"))

## plot
color <- c("slateblue", "violetred", "seagreen3", "steelblue1")
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- max(do2_scan) # overall maximum LOD score

for(j in 1:ceiling(ncol(do2_scan)/4)){
  if (j == ceiling(ncol(do2_scan)/4)) use = ((j-1)*4+1):ncol(do2_scan)
  else use = c(1:4) + 4*(j-1)
  
  for(i in 1:4){
    add=T
    if(i==1) add=F
    plot_scan1(do2_scan, gmap_list, lodcolumn=use[i], col=color[i], ylim=c(0, ymx*1.02), add=add)
  }
  
  legend("topleft", lwd=2, col=paste(color)[1:i], colnames(do2_scan)[use], bg="gray90")
  
}
pks <- find_peaks(do2_scan, gmap_list, threshold=4, drop=1.5)
hiPks <- pks[which(pks$lod > 30), ]
phenKeep <- unique(hiPks$lodcolumn)[1:4]
use <- which(colnames(do2_scan) %in% phenKeep)

for(i in 1:length(phenKeep)){
  add=T
  if(i==1) add=F
  plot_scan1(do2_scan, gmap_list, lodcolumn=use[i], col=color[i], ylim=c(0, ymx*1.02), add=add)
  legend("topleft", lwd=2, col=paste(color)[1:i], colnames(do2_scan)[use], bg="gray90")
}

