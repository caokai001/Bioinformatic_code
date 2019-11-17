- 文章：https://www.nature.com/articles/s41467-018-03828-2#Fig4
![](https://upload-images.jianshu.io/upload_images/9589088-6dd17a55d49591a5.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
- 具体介绍见星球:https://t.zsxq.com/zVVvVBI

```R
## Figure 4 - CBS and CTCF motifs


####################################
##R 安装包
####################################
source("http://www.bioconductor.org/biocLite.R")

if (!require('GenomicRanges')) biocLite(c("GenomicRanges"))
biocLite(c("rtracklayer", "AnnotationHub", "Rsamtools"))
if (!require('maftools')) biocLite(c("maftools"))

library(rtracklayer)


####################################
ctcf.motif=read.table("fimo_all.txt",sep="\t")
unique(ctcf.motif$V5-ctcf.motif$V4)+1 # 19

ctcf.motif=GRanges(seqnames=ctcf.motif$V3,
				IRanges(start=ctcf.motif$V4,end=ctcf.motif$V5),
				pval=ctcf.motif$V8,
				qval=ctcf.motif$V9,
				dir=ctcf.motif$V6,
				motif=ctcf.motif$V10)
length(ctcf.motif)   # [1] 192235


#####################################
## bedtools
#####################################
cat fimo_all.txt |awk 'BEGIN{FS=OFS="\t"}{print $3,$4,$5}' >fimo_all.bed

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome

bedtools slop -i fimo_all.bed -g hg19.genome -b 1000 >fimo_all_flank_1k.bed

## 修改snp格式转换成bed
cut -f 1-3 41467_2018_3828_MOESM6_ESM.txt|sed 1d >snp.pos.bed
## 取交集
intersectBed -a fimo_all_flank_1k.bed -b snp.pos.bed  -wa -wb|cut -f 4-6|sort -k1,1 -k2,2n >snp.pos.filter.bed

######
intersectBed -a fimo_all_flank_1k.bed -b snp.pos.bed  -wa -wb|sort -k1,1 -k2,2n  >test.txt

cat test.txt |awk 'BEGIN{FS=OFS="\t"}{print $0,$5-($2+1000),$5-($3-1000)}'|awk 'BEGIN{FS=OFS="\t"}{if ($7<0 && $8<0){print $0,$7}else if ($7>0 && $8 >0){print $0,$8}else{print $0,0}}' >snp_relative.bed

### 画图
cut -f 9 snp_relative.bed |sort |uniq -c|awk 'BEGIN{OFS="\t"}{print $2,$1}'|sort -k1,1n
 >plot_data.bed
 
>>>R 
library(tidyverse)
D=tbl_df(read.table("plot_data.bed",sep="\t"))

ggplot(data=D,aes(x=V1,y=V2))+geom_line()+ggtitle(paste("CTCF_context and somatic loci","chr=1,2,3",sep="\n"))+
	xlab("CBS + flanking region (bp)")+ylab("somatic count in each loci")+theme_classic()


#######
## 以motif 中心为0
#######

cat test.txt |awk '{FS=OFS="\t"}{print $0,$5-(($2+$3)/2)}'|cut -f 7|sort |uniq -c|awk 'BEGIN{OFS="\t"}{print $2,$1}'|sort -k1,1n >plot_data_notmerge0.bed

# 画图
D=tbl_df(read.table("plot_data_notmerge0.bed",sep="\t"))
p<-ggplot(data=D,aes(x=V1,y=V2))+geom_line()+ggtitle(paste("CTCF_context and somatic loci","chr=1,2,3",sep="\n"))+
+ xlab("CBS + flanking region (bp)")+ylab("somatic count in each loci")+theme_classic()


###########################################################################
# 下载Dnase roadmap  http://www.roadmapepigenomics.org/data/tables/all#
##       https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
# 下载CTCF
# https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=TF+ChIP-seq&award.project=ENCODE&assembly=hg19&target.investigated_as=transcription+factor&target.label=CTCF&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&advancedQuery=%40type%3AExperiment+date_released%3A%5B2009-01-01+TO+2018-04-30%5D&biosample_ontology.term_name=stomach
# https://www.encodeproject.org/experiments/ENCSR618QYE/  去除input
##########################################################################
## 下载基因组
cd /public/home/kcao/Desktop/new
tar -zxvf  /public/home/kcao/genome_human/hg19/chromFa.tar.gz
ls |grep "_" |xargs -i rm -rf {}
## fimo 拼接起开
cat output_chr1/fimo.tsv |grep -v "^$"|grep -v "#"  >c1.txt
# 其他染色体
ls|grep -v "fa"|grep -v -w "output_chr1"|sort -V|while read id;do echo "grep "MA" ${id}/fimo.tsv >c${id/output_chr/}.txt"|qsub -d ./ -N ${id} ;done
## for循环合并
for i in {1..22} X Y;do echo c$i.txt;done | xargs -i cat {} >> mynew.txt
[kcao@login new]$ sz mynew.txt
#############
# 取交集
## 传到工作站上面
R 
library("rtracklayer")
ctcf.motif=read.table("mynew.txt",sep="\t",header=1) 


ctcf.motif=GRanges(seqnames=ctcf.motif$sequence_name,
				IRanges(start=ctcf.motif$start,end=ctcf.motif$stop),
				pval=ctcf.motif$p.value,
				qval=ctcf.motif$q.value,
				dir=ctcf.motif$strand,
				motif=ctcf.motif$matched_sequence)
length(ctcf.motif)   # [1] 192235
unique(width(ctcf.motif)) # 19

##dnase.peak
dnase.peak=read.table("E094-DNase.all.peaks.v2.bed",sep="\t") 
dnase.peak=GRanges(seqnames=dnase.peak$V1,
				IRanges(start=dnase.peak$V2,end=dnase.peak$V3),
				score=dnase.peak$V5)

> sum(dnase.peak$score>0)
[1] 165200

## fdr0.01.peaks
dnase.fdr=read.table("E094-DNase.fdr0.01.peaks.v2.bed")
dnase.fdr=GRanges(seqnames=dnase.fdr$V1,IRanges(start=dnase.fdr$V2,end=dnase.fdr$V3))

## macs2.narrowpeak
dnase.mac=read.table("E094-DNase.macs2.narrowPeak")
dnase.mac=GRanges(seqnames=dnase.mac$V1,IRanges(start=dnase.mac$V2,dnase.mac$V3))

dnase=c(dnase.peak,dnase.fdr,dnase.mac)
dnase=reduce(dnase)
length(dnase)

> length(dnase)
[1] 223139

#### snptoGrange
somatic_gastric=read.table("41467_2018_3828_MOESM6_ESM.txt",header=T)

somatic_gastric<-GRanges(seqnames=somatic_gastric$seqnames,IRanges(start=somatic_gastric$start,somatic_gastric$end),
+ ral=somatic_gastric$ral,tal=somatic_gastric$tal,sid=somatic_gastric$sid,batch=somatic_gastric$batch)

somatic_gastric$id=c(1:length(somatic_gastric))
> length(somatic_gastric) 
[1] 4119812



z=findOverlaps(dnase,ctcf.motif)
dnase.ovl=dnase[queryHits(z)]
dnase.ovl=as.data.frame(dnase.ovl)
dnase.ovl=unique(dnase.ovl) 
motif.ovl=ctcf.motif[subjectHits(z)]
motif.ovl=as.data.frame(motif.ovl)
motif.ovl=unique(motif.ovl) 

# save ctcf union motifs to bed file
write.table(motif.ovl[,c("seqnames","start","end")],file="ctcf_motif_union.bed")

# count the number of somatic mutations at each site
site=GRanges(seqnames=motif.ovl$seqnames,IRanges(start=(motif.ovl$start+motif.ovl$end)/2,end=(motif.ovl$start+motif.ovl$end)/2),pval=motif.ovl$pval,dir=motif.ovl$dir,motif=motif.ovl$motif)

site=site+1000
unique(width(site))

> unique(width(site))
[1] 2001


> length(site)
[1] 84458

mut=findOverlaps(site,somatic_gastric) 
length(mut)

### --------
mut.table=somatic_gastric[subjectHits(mut)] 
mut.table$region.start=start(site[queryHits(mut)])
mut.table$region.end=end(site[queryHits(mut)])
mut.table$region.chr=seqnames(site[queryHits(mut)])
mut.table$pval=site[queryHits(mut)]$pval
mut.table$dir=site[queryHits(mut)]$dir
mut.table$motif=site[queryHits(mut)]$motif
length(unique(mut.table$id)) # 65966

mut.table$id=c(1:length(mut.table))

########################
mut.table=as.data.frame(mut.table)
mut.table=mut.table[order(mut.table$sid,mut.table$pval),]
mut.table=mut.table[!duplicated(mut.table$id),] 
mut.table$seqnames=as.character(mut.table$seqnames)
mut.table$region.chr=as.character(mut.table$region.chr)
sum(mut.table$seqnames==mut.table$region.chr) 


> sum(mut.table$seqnames==mut.table$region.chr) 
[1] 161097

#############
mut.table$motif.start=(mut.table$region.start+mut.table$region.end)/2 # 187 unique sids
overall.table=mut.table

# split mut.table by subtype
####-----------------------------
## 参考：https://blog.csdn.net/tb3039450/article/details/52557200
install.packages("readxl")
library(readxl)
subtype_classification <-read_excel("41467_2018_3828_MOESM4_ESM.xlsx",sheet=1,na="NA",skip=2)




GS=subtype_classification[which(subtype_classification$`Molecular Subtype`=="GS"),]$"Sample ID"
GS=c(GS,"apollo1_new") # 20


CIN=subtype_classification[which(subtype_classification$`Molecular Subtype`=="CIN"),]$"Sample ID" # 42
EBV=subtype_classification[which(subtype_classification$`Molecular Subtype`=="EBV"),]$"Sample ID" # 17
MSI=subtype_classification[which(subtype_classification$`Molecular Subtype`=="MSI"),]$"Sample ID" # 18

MSI=c(MSI,"CGP_donor_GC00031") # 19

######
## %in%  用于字符串包含查看
######

GS=GS[which(GS %in% unique(as.character(somatic_gastric$sid)))] # 11
CIN=CIN[which(CIN %in% unique(as.character(somatic_gastric$sid)))] # 41
EBV=EBV[which(EBV %in% unique(as.character(somatic_gastric$sid)))] # 17
MSI=MSI[which(MSI %in% unique(as.character(somatic_gastric$sid)))] # 19

> print(c(length(GS),length(CIN),length(EBV),length(MSI)))
[1] 11 41 17 20

> sum(unique(mut.table$sid) %in% c(MSI,EBV,GS,CIN))
[1] 88


mut.table=mut.table[which(mut.table$sid  %in% c(MSI,EBV,GS,CIN)),] 
# GS
gs=mut.table[which(mut.table$sid %in% GS),]
> nrow(gs)
[1] 3644

# MSI
msi=mut.table[which(mut.table$sid %in% MSI),] # 35152
nrow(msi)

# EBV
ebv=mut.table[which(mut.table$sid %in% EBV),] # 3976
nrow(ebv)

# CIN
cin=mut.table[which(mut.table$sid %in% CIN),] # 9222
nrow(cin)




















##########################################################Figure A
######################################################CIN subtype

> table(cin$dir) 

    -     + 
10928 10850 

cin.pos=cin[which(cin$dir=="+"),] 
cin.neg=cin[which(cin$dir=="-"),] 

cin.pos$pos=cin.pos$start-cin.pos$motif.start
sum(cin.pos$pos==0)






cin.pos.df=aggregate(seqnames~pos,cin.pos,length)
colnames(cin.pos.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
cin.pos.df=merge(df,cin.pos.df,all.x=TRUE,by="pos")
cin.pos.df$count=ifelse(is.na(cin.pos.df$count),0,cin.pos.df$count)

cin.neg$pos=cin.neg$motif.start-cin.neg$start
sum(cin.neg$pos==0) 

######
#画figA

cin.neg.df=aggregate(seqnames~pos,cin.neg,length)
colnames(cin.neg.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
cin.neg.df=merge(df,cin.neg.df,all.x=TRUE,by="pos")
cin.neg.df$count=ifelse(is.na(cin.neg.df$count),0,cin.neg.df$count)

cin.z=merge(cin.pos.df,cin.neg.df,by="pos")
cin.z$mut=cin.z$count.x+cin.z$count.y
cin.z$norm.mut=(cin.z$count.x+cin.z$count.y)/length(CIN)

p<-ggplot(data=cin.z,aes(x=pos,y=norm.mut))+
  geom_line()+
  ylim(c(0,2.5))+ggtitle(paste("CTCF_CIN","n=41",sep="\n"))+
  xlab("CBS + flanking region (bp)")+
  ylab("Normalized somatic substitutions")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line=element_line(colour="black"))

ggsave("CIN_without_ctcf_somatic_pos.pdf",p) 

######
#画figB

table(ebv$dir) 

> table(ebv$dir) 

   -    + 
4689 4607


ebv.pos=ebv[which(ebv$dir=="+"),] 
ebv.neg=ebv[which(ebv$dir=="-"),] 

ebv.pos$pos=ebv.pos$start-ebv.pos$motif.start
sum(ebv.pos$pos==0) 


ebv.pos.df=aggregate(seqnames~pos,ebv.pos,length)
colnames(ebv.pos.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
ebv.pos.df=merge(df,ebv.pos.df,all.x=TRUE,by="pos")
ebv.pos.df$count=ifelse(is.na(ebv.pos.df$count),0,ebv.pos.df$count)

ebv.neg$pos=ebv.neg$motif.start-ebv.neg$start
sum(ebv.neg$pos==0) 





ebv.neg.df=aggregate(seqnames~pos,ebv.neg,length)
colnames(ebv.neg.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
ebv.neg.df=merge(df,ebv.neg.df,all.x=TRUE,by="pos")
ebv.neg.df$count=ifelse(is.na(ebv.neg.df$count),0,ebv.neg.df$count)

ebv.z=merge(ebv.pos.df,ebv.neg.df,by="pos")
ebv.z$mut=ebv.z$count.x+ebv.z$count.y
ebv.z$norm.mut=(ebv.z$count.x+ebv.z$count.y)/length(EBV)

p2<-ggplot(data=ebv.z,aes(x=pos,y=norm.mut))+
  geom_line()+
  ylim(c(0,2.5))+
  ggtitle(paste("CTCF_EBV","n=17",sep="\n"))+
  xlab("CBS + flanking region (bp)")+
  ylab("Normalized somatic substitutions")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"))

ggsave("EBV_without_ctcf_somatic_pos.pdf",p2) 


### Fig C
## GS

> table(gs$dir)

   -    + 
1824 1820

gs.pos=gs[which(gs$dir=="+"),] 
gs.neg=gs[which(gs$dir=="-"),] 

gs.pos$pos=gs.pos$start-gs.pos$motif.start
sum(gs.pos$pos==0)


gs.pos.df=aggregate(seqnames~pos,gs.pos,length)
colnames(gs.pos.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
gs.pos.df=merge(df,gs.pos.df,all.x=TRUE,by="pos")
gs.pos.df$count=ifelse(is.na(gs.pos.df$count),0,gs.pos.df$count)

gs.neg$pos=gs.neg$motif.start-gs.neg$start
sum(gs.neg$pos==0)


gs.neg.df=aggregate(seqnames~pos,gs.neg,length)
colnames(gs.neg.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
gs.neg.df=merge(df,gs.neg.df,all.x=TRUE,by="pos")
gs.neg.df$count=ifelse(is.na(gs.neg.df$count),0,gs.neg.df$count)

gs.z=merge(gs.pos.df,gs.neg.df,by="pos")
gs.z$mut=gs.z$count.x+gs.z$count.y
gs.z$norm.mut=(gs.z$count.x+gs.z$count.y)/length(GS)

p3<-ggplot(data=gs.z,aes(x=pos,y=norm.mut))+
  geom_line()+
  ylim(c(0,2.5))+
  ggtitle(paste("CTCF_GS","n=11",sep="\n"))+
  xlab("CBS + flanking region (bp)")+
  ylab("Normalized somatic substitutions")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        axis.line=element_line(colour="black"))

ggsave("GS_without_ctcf_somatic_pos.pdf",p3) 




### Fig D 
## MSI
> table(msi$dir) 

    -     + 
44564 44415 


msi.pos=msi[which(msi$dir=="+"),] 
msi.neg=msi[which(msi$dir=="-"),]

msi.pos$pos=msi.pos$start-msi.pos$motif.start
sum(msi.pos$pos==0)




msi.pos.df=aggregate(seqnames~pos,msi.pos,length)
colnames(msi.pos.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
msi.pos.df=merge(df,msi.pos.df,all.x=TRUE,by="pos")
msi.pos.df$count=ifelse(is.na(msi.pos.df$count),0,msi.pos.df$count)

msi.neg$pos=msi.neg$motif.start-msi.neg$start
sum(msi.neg$pos==0)




msi.neg.df=aggregate(seqnames~pos,msi.neg,length)
colnames(msi.neg.df)[2]="count"
df=data.frame(pos=seq(-1000,1000,1),x=0)
msi.neg.df=merge(df,msi.neg.df,all.x=TRUE,by="pos")
msi.neg.df$count=ifelse(is.na(msi.neg.df$count),0,msi.neg.df$count)

msi.z=merge(msi.pos.df,msi.neg.df,by="pos")
msi.z$mut=msi.z$count.x+msi.z$count.y
msi.z$norm.mut=(msi.z$count.x+msi.z$count.y)/length(MSI)

P4<-ggplot(data=msi.z,aes(x=pos,y=norm.mut))+
  geom_line()+
  ylim(c(0,2.5))+
  ggtitle(paste("CTCF_MSI","n=19",sep="\n"))+
  xlab("CBS + flanking region (bp)")+
  ylab("Normalized somatic substitutions")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave("MSI_without_ctcf_somatic_pos.pdf",P4) 



















################
### fig F 
################
# CTCF motifs + flank 5bp
# motif.ovl表示过滤掉Dnase外部的motif[文章中过滤了ctcfpeak]
motifs=GRanges(seqnames=motif.ovl$seqnames,
               IRanges(motif.ovl$start,motif.ovl$end),
               dir=motif.ovl$dir,
               pval=motif.ovl$pval,
               motif=motif.ovl$motif)

> motifs=motifs+5
> length(motifs)
[1] 84458


z=findOverlaps(somatic_gastric,motifs)
t=as.data.frame(somatic_gastric[queryHits(z)]) # 3463
t$motif.seq=as.character(seqnames(motifs[subjectHits(z)]))
t$motif.start=start(motifs[subjectHits(z)])
t$motif.end=end(motifs[subjectHits(z)])
t$dir=motifs[subjectHits(z)]$dir
t$pval=motifs[subjectHits(z)]$pval
t$motif=motifs[subjectHits(z)]$motif



t=t[order(t$id,t$pval),]
t=t[!duplicated(t$id),] # 3093
t.pos=t[which(t$dir=="+"),]
t.neg=t[which(t$dir=="-"),]


################### FigF：Part2
## All CBS mutations Reference Alignment

# All
ref.seq.pos=GRanges(seqnames=t.pos$seqnames,
                IRanges(t.pos$motif.start,t.pos$motif.end),
                dir=t.pos$dir,
                motif=t.pos$motif,
                mut=t.pos$start,
                alt=t.pos$tal,
                ref=t.pos$ral)
				
				
				
library(BSgenome.Hsapiens.UCSC.hg19)
extrSeq=Views(Hsapiens,ref.seq.pos)

z=as(extrSeq,"DNAStringSet") 
dnaSeq.pos=as.data.frame(ref.seq.pos)
dnaSeq.pos$dna=as.character(z) # 2001


dnaSeq.pos$motif=toupper(dnaSeq.pos$motif)
sum(dnaSeq.pos$motif==dnaSeq.pos$dna)






#####################
# 画图
#####################

dnaSeq.pos$a1=substr(dnaSeq.pos$dna,1,1)
dnaSeq.pos$a2=substr(dnaSeq.pos$dna,2,2)
dnaSeq.pos$a3=substr(dnaSeq.pos$dna,3,3)
dnaSeq.pos$a4=substr(dnaSeq.pos$dna,4,4)
dnaSeq.pos$a5=substr(dnaSeq.pos$dna,5,5)
dnaSeq.pos$a6=substr(dnaSeq.pos$dna,6,6)
dnaSeq.pos$a7=substr(dnaSeq.pos$dna,7,7)
dnaSeq.pos$a8=substr(dnaSeq.pos$dna,8,8)
dnaSeq.pos$a9=substr(dnaSeq.pos$dna,9,9)
dnaSeq.pos$a10=substr(dnaSeq.pos$dna,10,10)
dnaSeq.pos$a11=substr(dnaSeq.pos$dna,11,11)
dnaSeq.pos$a12=substr(dnaSeq.pos$dna,12,12)
dnaSeq.pos$a13=substr(dnaSeq.pos$dna,13,13)
dnaSeq.pos$a14=substr(dnaSeq.pos$dna,14,14)
dnaSeq.pos$a15=substr(dnaSeq.pos$dna,15,15)
dnaSeq.pos$a16=substr(dnaSeq.pos$dna,16,16)
dnaSeq.pos$a17=substr(dnaSeq.pos$dna,17,17)
dnaSeq.pos$a18=substr(dnaSeq.pos$dna,18,18)
dnaSeq.pos$a19=substr(dnaSeq.pos$dna,19,19)
dnaSeq.pos$a20=substr(dnaSeq.pos$dna,20,20)
dnaSeq.pos$a21=substr(dnaSeq.pos$dna,21,21)
dnaSeq.pos$a22=substr(dnaSeq.pos$dna,22,22)
dnaSeq.pos$a23=substr(dnaSeq.pos$dna,23,23)
dnaSeq.pos$a24=substr(dnaSeq.pos$dna,24,24)
dnaSeq.pos$a25=substr(dnaSeq.pos$dna,25,25)
dnaSeq.pos$a26=substr(dnaSeq.pos$dna,26,26)
dnaSeq.pos$a27=substr(dnaSeq.pos$dna,27,27)
dnaSeq.pos$a28=substr(dnaSeq.pos$dna,28,28)
dnaSeq.pos$a29=substr(dnaSeq.pos$dna,29,29)





####################
## t.neg
####################
ref.seq.neg=GRanges(seqnames=t.neg$seqnames,
                    IRanges(t.neg$motif.start,t.neg$motif.end),
                    dir=t.neg$dir,
                    motif=t.neg$motif,
                    mut=t.neg$start,
                    alt=t.neg$tal,
                    ref=t.neg$ral)
extrSeq=Views(Hsapiens,ref.seq.neg)

z=as(extrSeq,"DNAStringSet") 
dnaSeq.neg=as.data.frame(ref.seq.neg)
dnaSeq.neg$dna=as.character(z) # 1542

dnaSeq.neg$motif=toupper(dnaSeq.neg$motif)
sum(dnaSeq.neg$motif==dnaSeq.neg$dna) # 0




z=reverseComplement(z)
dnaSeq.neg$rev.dna=as.character(z)
sum(dnaSeq.neg$motif==dnaSeq.neg$rev.dna)


dnaSeq.neg$a1=substr(dnaSeq.neg$rev.dna,1,1)
dnaSeq.neg$a2=substr(dnaSeq.neg$rev.dna,2,2)
dnaSeq.neg$a3=substr(dnaSeq.neg$rev.dna,3,3)
dnaSeq.neg$a4=substr(dnaSeq.neg$rev.dna,4,4)
dnaSeq.neg$a5=substr(dnaSeq.neg$rev.dna,5,5)
dnaSeq.neg$a6=substr(dnaSeq.neg$rev.dna,6,6)
dnaSeq.neg$a7=substr(dnaSeq.neg$rev.dna,7,7)
dnaSeq.neg$a8=substr(dnaSeq.neg$rev.dna,8,8)
dnaSeq.neg$a9=substr(dnaSeq.neg$rev.dna,9,9)
dnaSeq.neg$a10=substr(dnaSeq.neg$rev.dna,10,10)
dnaSeq.neg$a11=substr(dnaSeq.neg$rev.dna,11,11)
dnaSeq.neg$a12=substr(dnaSeq.neg$rev.dna,12,12)
dnaSeq.neg$a13=substr(dnaSeq.neg$rev.dna,13,13)
dnaSeq.neg$a14=substr(dnaSeq.neg$rev.dna,14,14)
dnaSeq.neg$a15=substr(dnaSeq.neg$rev.dna,15,15)
dnaSeq.neg$a16=substr(dnaSeq.neg$rev.dna,16,16)
dnaSeq.neg$a17=substr(dnaSeq.neg$rev.dna,17,17)
dnaSeq.neg$a18=substr(dnaSeq.neg$rev.dna,18,18)
dnaSeq.neg$a19=substr(dnaSeq.neg$rev.dna,19,19)
dnaSeq.neg$a20=substr(dnaSeq.neg$rev.dna,20,20)
dnaSeq.neg$a21=substr(dnaSeq.neg$rev.dna,21,21)
dnaSeq.neg$a22=substr(dnaSeq.neg$rev.dna,22,22)
dnaSeq.neg$a23=substr(dnaSeq.neg$rev.dna,23,23)
dnaSeq.neg$a24=substr(dnaSeq.neg$rev.dna,24,24)
dnaSeq.neg$a25=substr(dnaSeq.neg$rev.dna,25,25)
dnaSeq.neg$a26=substr(dnaSeq.neg$rev.dna,26,26)
dnaSeq.neg$a27=substr(dnaSeq.neg$rev.dna,27,27)
dnaSeq.neg$a28=substr(dnaSeq.neg$rev.dna,28,28)
dnaSeq.neg$a29=substr(dnaSeq.neg$rev.dna,29,29)


> df2=rbind(dnaSeq.pos[,12:40],dnaSeq.neg[,13:41])
> dim(df2)
[1] 3984   29

df=matrix(0,nrow=4,ncol=29)

rownames(df)=c("A","C","G","T")
for (i in 1:29){
  print(i)
  df[names(table(df2[,i])),i]=as.numeric(table(df2[,i]))
}


colnames(df)=colnames(df2)
print(df)

BiocManager::install("seqLogo")
library(seqLogo)
pdf()
pwm=makePWM(as.matrix(df/unique(colSums(df))),alphabet="DNA")

> pdf("CTCF_gs_seqlogo.pdf")
> seqLogo(pwm,ic.scale=TRUE,xaxis=TRUE,yaxis=TRUE,xfontsize=15,yfontsize=15)
> dev.off()




#####
## Alternate Alignment  :
## 突变位置替换后序列seqlogo 图
#####
dnaSeq.mut.pos=dnaSeq.pos
dnaSeq.mut.pos$alt=as.character(dnaSeq.mut.pos$alt)
dnaSeq.mut.pos$num=dnaSeq.mut.pos$mut-dnaSeq.mut.pos$start+1
for (i in 1:nrow(dnaSeq.mut.pos)){
  dnaSeq.mut.pos[i,paste("a",dnaSeq.mut.pos[i,"num"],sep="")]<-dnaSeq.mut.pos[i,"alt"]
}

dnaSeq.mut.neg=dnaSeq.neg
dnaSeq.mut.neg$alt=as.character(dnaSeq.mut.neg$alt)
dnaSeq.mut.neg$num=dnaSeq.mut.neg$end-dnaSeq.mut.neg$mut+1
for (i in 1:nrow(dnaSeq.mut.neg)){
  dnaSeq.mut.neg[i,paste("a",dnaSeq.mut.neg[i,"num"],sep="")]<-as.character(reverseComplement(DNAString(dnaSeq.mut.neg[i,"alt"])))
}



df7=rbind(dnaSeq.mut.pos[,12:40],dnaSeq.mut.neg[,13:41]) # 3093
df=matrix(0,nrow=4,ncol=29)
rownames(df)=c("A","C","G","T")
for (i in 1:29){
  print(i)
  df[names(table(df7[,i])),i]=as.numeric(table(df7[,i]))
}


colnames(df)=colnames(df7)
print(df)

pwm=makePWM(as.matrix(df/unique(colSums(df))),alphabet="DNA")
pdf("CTCF_mut_gs_seqlogo.pdf")
seqLogo(pwm,ic.scale=TRUE,xaxis=TRUE,yaxis=TRUE,xfontsize=15,yfontsize=15)
dev.off()





df=rbind(dnaSeq.mut.pos[,c("ref","alt","num")],dnaSeq.mut.neg[,c("ref","alt","num")]) # 3093
df$count=1
df=aggregate(count~ref+alt+num,df,sum) # 317

for (i in 1:nrow(df)){
  if (!df$ref[i] %in% c("C","T")){
    df$alt[i]=as.character(reverseComplement(DNAString(df$alt[i])))
    df$ref[i]=as.character(reverseComplement(DNAString(df$ref[i])))
  }
}
df$mut=paste(df$ref,df$alt,sep="")
df$mut=factor(df$mut,levels=c("TG","TC","TA","CT","CG","CA"))
df$num=factor(df$num,levels=c(1:29))

df=aggregate(count~ref+alt+num+mut,df,sum) # 169

p5<-ggplot(df,aes(x=num,y=count,fill=mut))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = c("#33CC00","#33CCFF","#FF9933","#CC33FF","#FFFF33","#FF0000"))

ggsave("mutate_count_motif_flank_5bp.pdf",p5)

```



![结果汇报](https://upload-images.jianshu.io/upload_images/9589088-6cbed9895118a3ee.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
