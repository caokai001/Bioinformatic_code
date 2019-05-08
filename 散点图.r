setwd("C:/Users/16926/Desktop/111/")
par(mar = c(5, 5, 5, 3),cex.lab  = 2,mfcol=c(3,4))
type<-c("27AC-siEAF1","27AC-siEAF2","27AC-siGFP",
        "k27me3-siEAF1","k27me3-siEAF2","k27me3-siGFP",
        "k4me1-siEAF1","k4me1-siEAF2","k4me1-siGFP",
        "K4me3-siEAF1","K4me3-siEAF2","K4me3-siGFP"
        )
for(i in c(1:12)){
  data<-read.table(paste(type[i],".tab",sep = ""),header = T)
  cols <- densCols(data[,c(4,5)], colramp=colorRampPalette(
    c("DarkBlue","blue","cyan","green","yellow","orange", "red")), nbin = 100)
  plot(data[,c(4,5)], 
       xlab="replicate_1",ylab="replicate_2",main=type[i],
       col=cols,pch=20,cex.main=2.2)
  legend("topleft", cex=2,bty="n",
         paste("Pcc = ", round(cor(data[,c(4,5)], method = "pearson")[1,2],3)))
}
