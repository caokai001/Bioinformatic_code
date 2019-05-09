######################deeptools 里面获得tab 文件，10k 为bin,画散点图。
library(stringi)
file=read.table("si_unique_bam_readCounts.tab",header = T)

A=names(file)
#numbers=c(10,12,14,20,22,24,26)
numbers=c(4,6,8,10,12,14,16,18,20,22,24,26)
for (i in numbers){
  png(paste0(stri_sub(A[i],,-15),"_correlation.png"),width = 480, height = 480)
  data<-file[,i:(i+1)]   ###add ()
  d<-log2(data)  ###+1
  cols <- densCols(d[,c(1,2)], colramp=colorRampPalette(
    c("DarkBlue","blue","cyan","green","yellow","orange", "red")), nbin = 100)
  plot(d[,c(1,2)], 
       xlab=A[i],ylab=A[i+1],main=stri_sub(A[i],,-15),
       col=cols,pch=20,cex.main=1.5,cex.lab=1.3)
  legend("topleft", cex=1,bty="n",
         paste("Spearman = ", round(cor(data[,c(1,2)], method = "spearman")[1,2],3)))
  legend("topright", cex=1,bty="n",
         paste("pcc = ", round(cor(data[,c(1,2)], method = "pearson")[1,2],3)))
  dev.off()
}

###效果如下：https://github.com/caokai001/Bioinformatic_code/blob/master/%E5%9B%BE%E7%89%87%E8%B5%84%E6%96%99%E5%BA%93/%E6%95%A3%E7%82%B9%E5%9B%BE.png
