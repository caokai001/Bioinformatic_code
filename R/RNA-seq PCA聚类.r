myfpkm<-read.table("All_gene_fpkm.xls",header=TRUE,comment.char="",sep = "\t",check.names=FALSE,row.names=1)  ##基因名为行名
probesetvar = apply(myfpkm, 1, var)
ord = order(probesetvar, decreasing=TRUE)[1:200]       ##ord 函数确定排序后的数值原来的行数，index
pca = prcomp(t(myfpkm[ord,]), scale=TRUE)
ss=summary(pca)

#绘图：

plot(pca$x[,1:2],col=rep(c(1,2,3,4,1,2,3,4),each=3),pch=rep(c(16,17),each=12))


#或者3D：

library(scatterplot3d)
scatterplot3d(pca$x[,1:3],color=rep(c(1,2,3,4,1,2,3,4),each=3),pch=rep(c(16,17),each=12))
