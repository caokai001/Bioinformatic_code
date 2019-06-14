## 1.看CTCF中snp 突变位置，与peak 关系。已知peak区域motif ,再和snp 关联，大部分都在非关键位置突变
### motif 上snp 个数
```
snp_pos=read.table("LNCaP_snp_count")
motif=substring("TGGCCACCAGGGGGCGCTA", 1:19, 1:19)
colnames(snp_pos)=c("count","pos")
snp_pos$pos=as.factor(snp_pos$pos)
ggplot(snp_pos,aes(x=pos,y=count))+geom_bar(stat="identity")+
  ggtitle("LNCaP CTCF motif and snp")+scale_x_discrete(labels=motif)
```
![pos 与snp](https://upload-images.jianshu.io/upload_images/9589088-5b480f9ae7a62409.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
###
```
LNCaP_judge=read.table("LNCaP_judge.txt")
motif=substring("TGGCCACCAGGGGGCGCTA", 1:19, 1:19)
tmp=LNCaP_judge[,12:15]
colnames(tmp)=c("ref","mut","pos","judge")
tmp$pos=as.factor(tmp$pos)
ggplot(tmp,aes(x=pos,fill=judge))+geom_bar()+ggtitle("LNCaP CTCF motif and snp")+
  labs(fill="Ref base")+scale_x_discrete(labels=motif)+xlab("motif pos")+coord_flip()
```
![yes 突变成motif 一样序列](https://upload-images.jianshu.io/upload_images/9589088-d2fbd73a9bb9423a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

>统计bar图：
1.需要考虑（-）情况
2.杂合突变
3.一个peak对应多个motif

## 修改
