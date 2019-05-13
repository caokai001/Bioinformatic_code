###练习 注释及富集    https://zhuanlan.zhihu.com/p/35510434
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(VennDiagram)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
setwd("C:\\Users\\16926\\Desktop\\第三次-bioclass")
# LNCaP_CTRL_R1881<-readPeakFile("ChIP-seq_LNCaP_CTRL_R1881.narrowPeak")
# LNCaP_EZH2_HA_rep1<-readPeakFile("ChIP-seq_LNCaP_EZH2-HA_rep1.narrowPeak")
# LNCaP_EZH2_HA_rep2<-readPeakFile("ChIP-seq_LNCaP_EZH2-HA_rep2.narrowPeak")
PLKO_EZH2<-readPeakFile("SRR6355936.uniq.bam.narrow_peaks.narrowPeak_1.bed")
PLKO_H3K27ac<-readPeakFile("SRR6355938.uniq.bam.narrow_peaks.narrowPeak_1.bed")
PLKO_H3K27me3<-readPeakFile("SRR6355940.uniq.bam.broad_peaks.broadPeak_1.bed")
###接着我们需要用peakAnno函数对其进行peak注释
PLKO_EZH2_peakAnno <- annotatePeak(PLKO_EZH2, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb="org.Hs.eg.db")
PLKO_H3K27ac_peakAnno <- annotatePeak(PLKO_H3K27ac, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb="org.Hs.eg.db")
PLKO_H3K27me3_peakAnno <- annotatePeak(PLKO_H3K27me3, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb="org.Hs.eg.db")
###注释到的gene，绘制venn图
venn_list <- list(EZH2=as.data.frame(PLKO_EZH2_peakAnno )$geneId,
                  H3K27ac=as.data.frame(PLKO_H3K27ac_peakAnno)$geneId,
                  H3K27me3=as.data.frame(PLKO_H3K27me3_peakAnno)$geneId)
venn.diagram(venn_list, imagetype = "png", fill = c("blue", "green","red"), 
             alpha = c(0.5, 0.5, 0.5)  ,filename = "venn.png")
#vennplot(venn_list)


###peak 注释可视化
covplot(PLKO_EZH2,weightCol = "V5")
covplot(PLKO_EZH2, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

###获取EZGH2附近的基因,,从而得到tagMatrix并绘图
EZH2_gene<-as.data.frame(PLKO_EZH2_peakAnno)$geneId

EZH2_df<-as.data.frame(PLKO_EZH2_peakAnno)
EZH2_peak<-GRanges(EZH2_df[,1:12])

H3K27ac_df<-as.data.frame(PLKO_H3K27ac_peakAnno)
H3K27ac_peak<-GRanges(H3K27ac_df[H3K27ac_df$geneId %in% EZH2_gene, 1:12])

H3K27me3_df<-as.data.frame(PLKO_H3K27me3_peakAnno)
H3K27me3_peak<-GRanges(H3K27me3_df[H3K27me3_df$geneId %in% EZH2_gene, 1:12])
###生成EZH2注释gene,出三种ChIP-seq的list文件
list_peak <- list(H3K27ac=H3K27ac_peak, H3K27ac=H3K27ac_peak, EZH2=EZH2_peak)


promoter <- promoters(txdb, upstream = 3000, downstream = 3000)
tagMatrix <- lapply(list_peak, getTagMatrix, window=promoter)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency")

tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color=NULL)








### ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrix <- getTagMatrix(PLKO_EZH2, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

##################################################################GOrich
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

keytypes(org.Hs.eg.db)
data<-read.table("EZH2_gene.list")
data$V1<-as.character(data$V1)
test1 = bitr(data$V1, fromType="ENSEMBLTRANS", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
head(test1,2)

###三、 GO分析
#3.1 groupGO 富集分析
ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",level = 3,readable = TRUE)

##3.2 enrichGO 富集分析
ego_BP <- enrichGO(gene = test1$ENTREZID, 
                #背景基因集
                OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
				        #keyType = 'ENSEMBL',
                ont = "BP", #也可以是 CC  BP  MF中的一种
                pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                qvalueCutoff = 1,
				        readable = TRUE)#Gene ID 转成gene Symbol ，易读



dotplot(ego_BP,title="EnrichmentGO_BP_dot")
barplot(ego_BP, showCategory=20,title="EnrichmentGO_BP")#条状图，按p从小到大排，绘制前20个Term
plotGOgraph(ego_BP)

#################################################dotplot
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
y=ego_BP@result
####
## 分别后去分号前面和后面的数，并变成数值
forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))
## 增加数值表示的一列GeneRatio
y$GeneRatio = forward/backward
showCategory =10
font.size =10

y %>% 
  ## 安装p值排序，选区既定数目的行
  arrange(p.adjust) %>% 
  slice(1:showCategory) %>% 
  ## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  #margin 上右下左
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=45))


#############################################barplot
## 设定显示的数目
showCategory =8
## 设定字体的大小
font.size =12
y %>% 
  ## 安装p值排序，选区既定数目的行
  arrange(p.adjust) %>% 
  slice(1:showCategory) %>% 
  ## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
  ggplot(aes(x=forcats::fct_reorder(Description,p.adjust,.desc = T),y=Count,fill=p.adjust))+ 
  ## 画出bar图
  geom_bar(stat="identity")+
  coord_flip()+
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 如果用ylab("")或出现左侧空白
  labs(x=NULL,y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))


