###用shell 写个循环跑起来
###usage : Rscript 脚本.r input.file output.file
library("gProfileR")
library("dplyr")
library("parallel")  ###并行计算
library ("plyr")
args=commandArgs(T)
input.file<-args[1]
output.file<-args[2]

A=read.table(input.file)
A=head(A$V1,20)
###定义转换ID函数
fun <- function(x){
  T<-gconvert(x, organism = "hsapiens", target="ENTREZGENE", filter_na =F)
  T<-tbl_df(T)#%>%arrange(alias.number)
  return(T)
}
##可用核心数目
cl.cores <- detectCores()
#利用并行计算
cl <- makeCluster(3)  # 初始化四核心集群
clusterEvalQ(cl,library("gProfileR"))
clusterEvalQ(cl,library("dplyr"))
results <- parLapply(cl,A,fun)#调用parLapply并行计算平方函数
df <- ldply (results, data.frame)#整合结果
stopCluster(cl) # 关闭集群
write.table(df,output.file,quote=FALSE,row.names=FALSE,sep=",")
