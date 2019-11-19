```
pwm <- function(freq, total, bg=0.25){
  #using the formulae above
  p <- (freq + (sqrt(total) * 1/4)) / (total + (4 * (sqrt(total) * 1/4)))
  log2(p/bg)
}


Cal_score<-function(seq,mm){
  seq<-toupper(seq)
  x <- strsplit(x=seq,split='')
  seq_score <- vector()
  for (i in 1:(dim(mm)[2]) ){
    seq_score[i] <- as.numeric(mm[x[[1]][i],i])
  }
  sum(seq_score)
}




## 1.1  拓宽PWM矩阵

df2=rbind(dnaSeq.pos[,12:40],dnaSeq.neg[,13:41]) # 3093
df=matrix(0,nrow=4,ncol=29)
rownames(df)=c("A","C","G","T")
for (i in 1:29){
  print(i)
  df[names(table(df2[,i])),i]=as.numeric(table(df2[,i]))
}

mm<-pwm(df,unique(colSums(df)))

#seq <- 'ACAGATGACCACTAGATGCCGCTGCCGGG'
Cal_score(seq,mm)





## 1.2 dnaSeq.pos

head(dnaSeq.pos)
dnaSeq.pos.score_test=dnaSeq.pos
head(dnaSeq.pos.score_test)

for (i in 1:nrow(dnaSeq.pos.score_test)){
dnaSeq.pos.score_test[i,"ref_score"]=Cal_score(dnaSeq.pos.score_test[i,"dna"],mm)
mut_seq=as.character(dnaSeq.pos.score_test[i,"dna"])
mut.relative.pos=dnaSeq.pos.score_test[i,"mut"]-dnaSeq.pos.score_test[i,"start"]+1
substr(mut_seq,mut.relative.pos,mut.relative.pos)<-as.character(dnaSeq.pos.score_test[,"alt"])[i]

dnaSeq.pos.score_test[i,"alt_score"]=Cal_score(mut_seq,mm) 
}

dnaSeq.pos.score_test$pos=dnaSeq.pos.score_test$mut-dnaSeq.pos.score_test$start +1






## 1.3 DNAseq.neg

head(dnaSeq.neg)
dnaSeq.neg.score_test=dnaSeq.neg
head(dnaSeq.neg.score_test)
i=1
for (i in 1:nrow(dnaSeq.neg.score_test)){
dnaSeq.neg.score_test[i,"ref_score"]=Cal_score(dnaSeq.neg.score[i,"rev.dna"],mm)
mut_seq=as.character(dnaSeq.neg.score_test[i,"rev.dna"])
mut.relative.pos=dnaSeq.neg.score_test[i,"mut"]-dnaSeq.neg.score_test[i,"start"]+1
substr(mut_seq,mut.relative.pos,mut.relative.pos)<-as.character(dnaSeq.neg.score_test[,"alt"])[i]

dnaSeq.neg.score_test[i,"alt_score"]=Cal_score(mut_seq,mm) 
}

dnaSeq.neg.score_test$pos=dnaSeq.neg.score_test$mut-dnaSeq.neg.score_test$start +1







## 1.4 画出

tmp_1=dnaSeq.pos.score_test[,c("ref","alt","ref_score","alt_score","pos","dir")]
tmp_2=dnaSeq.neg.score_test[,c("ref","alt","ref_score","alt_score","pos","dir")]
merge_score=rbind(tmp_1,tmp_2)



file.create("./ctcf_binding_affinity/wilcox_test.txt")
library(ggplot2)
library(ggpubr)

for (i in c(2,3,4,5,6,7,14)){
labels=c(-4,-3,-2,-1,1,2,9)
label_mut=as.character(labels[which(i==c(2,3,4,5,6,7,14))])


## R 将结果文件出到  https://www.cnblogs.com/timeisb

d_tmp=merge_score[merge_score$pos==i,]
d=gather(d_tmp,Type,value,-c(ref,alt,dir,pos))

my_comparisons <- list(c("ref_score", "alt_score"))
p2<-ggboxplot(d, x="Type", y="value", color = "Type", palette = c("#00AFBB", "#E7B800")),                             
add = "jitter", shape="Type")+
stat_compare_means(comparisons = my_comparisons,paired = TRUE,method = "wilcox.test") +
ggtitle(paste0("Score changed before and after the mutation\n mut_pos:",label_mut,"  n: ",nrow(d)))                  
#stat_compare_means(label.y = 20)   # 调整label.y 轴高度


ggsave(filename =paste0("ctcf_binding_affinity/mut_",i,".pdf"),p2,device="pdf")
#pdf(paste0("ctcf_binding_affinity/mut_",i,".pdf")) # , width=8, height=5
#p2
#dev.off()


## 配对样本t测验： https://zhuanlan.zhihu.com/p/36038655

## 单变量，多变量正态性检测 https://zhuanlan.zhihu.com/p/41110516

## R 将结果文件出到  https://www.cnblogs.com/timeisbiggestboss/p/7766115.html

#d_test<-d_tmp$ref_score-d_tmp$alt_score
#qqnorm(d_test);qqline(d_test)


sink("./ctcf_binding_affinity/wilcox_test.txt",append=TRUE)  # 输出重定向
print(paste0("#######################突变位置在",i,"处################################"))                        
print("####正态性")
print(shapiro.test(d_tmp$ref_score))
#shapiro.test(d_tmp$alt_score)
print("####配对的非参数检验")
print(wilcox.test(d_tmp$ref_score,d_tmp$alt_score,paired=TRUE))
sink()  # 结束重定向
}



```

>结果
![](https://upload-images.jianshu.io/upload_images/9589088-3499b9f8fe7c6a16.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
