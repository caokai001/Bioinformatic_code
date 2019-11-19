install.packages("grid")
install.packages("gridExtra")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("ppbeeswarm")  ##峰巢图
install.packages("ggbeeswarm")
install.packages("cowplot")
install.packages("pheatmap")
install.packages("ggsci")
install.packages("ggpubr")
library(grid)
library(gridExtra)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(ggpubr)


Sys.setlocale(category="LC_ALL",locale="en_US.UTF-8")
par(family='STXihei')
Sys.setlocale(locale="chinese")


#####1.JC 封面
#tips:Geom_text（）将文本直接添加到绘图中。 geom_label（）在文本后面绘制一个矩形，使其更易于阅读。
# 参考：https://www.plob.org/article/7735.html
## MY_code
pdf("JC_Bioinformatics_R_dxiong.pdf",family="GB1", width=8, height=5)
ggplot(data=data.frame())+
  geom_text(aes(x=1,y=1,label="JC Bioinformatics\nR语言 第二讲\nggplot2"),size=10)+
  geom_text(aes(x=1.5,y=0.2,label="信息学院\n2019.07.15\n"),size=5)+
  xlim(0,2)+ylim(0,1.5)+
  theme_bw()+
  labs(title = element_blank(),x=element_blank(),y=element_blank())+
  theme(panel.border = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank(),axis.text=element_blank())
dev.off()

################introduction  fig1###########
(ggplot(data=data.frame())
  +geom_text(aes(x=1,y=1,label="JC Bioinformatics\nR语言 第二讲\nggplot2"),size=10)
  +geom_text(aes(x=1.5,y=0.2,label="信息学院\n2019.07.15\n"),size=5)
  +xlim(0,2)+ylim(0,1.5)
  +theme_bw()
  +theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
  +theme(panel.border = element_blank(),panel.grid = element_blank())
)





####2.主要内容
# hjust和vjust的值限定在0-1之间。 0,0左底对齐
# hjust控制水平横轴，0表示右适应; 1表示左适应。
# vjust控制vertical纵轴，0表示底部适应; 1表示头部适应。 http://blog.leanote.com/post/shenweiyan/r-hjust-vjust
##My_code
pdf("2.test.pdf",family="GB1")
ggplot(data=data.frame())+
  xlim(0,1)+ylim(0,1)+
  geom_text(aes(x=0,y=1,label="主要内容"),size=10,hjust=0,vjust=1)+
  geom_text(aes(x=0,y=0.9,label="\n1、基础绘图\n    密度图、直方图、散点图、柱状图、饼图、误差棒图、\n    箱线图、小提琴图、蜜蜂图、热图\n2、分面\n3、图形组合\n4、样式精修\n5、添加显著性检验\n6、模式图绘制展示"),size=5,hjust=0,vjust=1)+
  theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())+
  theme(panel.border = element_blank(),panel.grid = element_blank())
dev.off()

################main text##############
(ggplot(data=data.frame())
 +geom_text(aes(x=0,y=1,label="主要内容"),size=10,hjust=0,vjust=1)
 +geom_text(aes(x=0,y=0.9,label="\n1、基础绘图\n    密度图、直方图、散点图、柱状图、饼图、误差棒图、\n    箱线图、小提琴图、蜜蜂图、热图\n2、分面\n3、图形组合\n4、样式精修\n5、添加显著性检验\n6、模式图绘制展示"),size=5,hjust=0,vjust=1)
 +xlim(0,1)+ylim(0,1)
 +theme_bw()
 +theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
 +theme(panel.border = element_blank(),panel.grid = element_blank())
)





###3.主要图形展示
# read.table() 默认没有header,默认\t
# 修改图形颜色  http://lilibei.net/2017/04/10/%E5%A6%82%E4%BD%95%E4%BD%BF%E7%94%A8ggplot2%E4%BF%AE%E6%94%B9%E9%A2%9C%E8%89%B2/
# http://www.rpubs.com/lihaoyi/156592
# My_code
background<-(theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5)))
fpkm<-read.table("data_for_JCB.txt",header=T,row.names = "gene")
fpkm2<-as.data.frame(t(fpkm))
# 上方四个图
P1<-ggplot(fpkm2)+geom_density(aes(Zm00001d013156),fill="grey",alpha=0.7)+ggtitle("geom_density")+background
P2<-ggplot(fpkm2)+geom_histogram(aes(Zm00001d006775),bins = 30)+labs(title = "geom_histogram")+background+scale_x_log10()
fpkm2$tissue<-row.names(fpkm2)
P3<-ggplot(fpkm2)+geom_point(aes(x=tissue,y=Zm00001d006775,colour=tissue))+labs(title = "geom_point")+scale_color_discrete(1:length(fpkm2$tissue))+theme(legend.position = "none")+background
P4<-ggplot(fpkm2)+geom_bar(aes(tissue,Zm00001d006775),stat="identity")+labs(title = "geom_bar")+background

#  下方四个图
P5<-ggplot(fpkm2)+geom_boxplot(aes(x="1",y=Zm00001d013156,fill="blue"))+geom_boxplot(aes("2",Zm00001d006775),fill="red")+labs(title = "geom_boxplot")+background
P6<-ggplot(fpkm2)+geom_violin(aes("Zm00001d013156",Zm00001d013156),fill="blue")+geom_violin(aes("Zm00001d006775",Zm00001d006775),fill="red")+labs(title = "geom_violin")+background


P7<-ggplot(fpkm2)+geom_beeswarm(aes("Zm00001d013156",Zm00001d013156,color="blue"))+geom_beeswarm(aes("Zm00001d006775",Zm00001d006775),color="blue")+labs(title = "geom_beeswarm")+background

# 变形reshape
data<-melt(fpkm2[c(1:10),c(1:10,401)],id.vars = "tissue")
names(data)<-c("tissue","gene","FPKM")
P8<-ggplot(data)+geom_tile(aes(gene,tissue,fill=log2(FPKM+1)))+labs(title = "geom_tile")+scale_fill_gradient2(high="red",low="blue",mid="white",midpoint = 6)+background

#grid.arrange(P1,P2,P3,P4,P5,P6,P7,P8,nrow=2)
plot_grid(P1,P2,P3,P4,P5,P6,P7,P8,labels = "AUTO",nrow = 2,hjust = 0.2)


################图形展示###################
background<-(theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none"))
fpkm<-read.table("data_for_JCB.txt",header=T,sep = "\t",row.names = "gene")
fpkm2<-as.data.frame(t(fpkm))
density_plot<-ggplot(fpkm2)+geom_density(aes(Zm00001d013156),fill="grey",alpha=0.7)+labs(title = "geom_density")+background
hist_plot<-ggplot(fpkm2)+geom_histogram(aes(Zm00001d006775),bins = 30)+labs(title = "geom_histogram")+background
fpkm2$tissue<-row.names(fpkm2)
point_plot<-ggplot(fpkm2)+geom_point(aes(tissue,Zm00001d006775))+labs(title = "geom_point")+background
bar_plot<-ggplot(fpkm2)+geom_bar(aes(tissue,Zm00001d006775),stat="identity")+labs(title = "geom_bar")+background
box_plot<-ggplot(fpkm2)+geom_boxplot(aes("Zm00001d013156",Zm00001d013156),fill="grey50")+geom_boxplot(aes("Zm00001d006775",Zm00001d006775),fill="grey50")+labs(title = "geom_boxplot")+background
violin_plot<-ggplot(fpkm2)+geom_violin(aes("Zm00001d013156",Zm00001d013156),fill="grey50")+geom_violin(aes("Zm00001d006775",Zm00001d006775),fill="grey50")+labs(title = "geom_violin")+background
bee_plot<-ggplot(fpkm2)+geom_beeswarm(aes("Zm00001d013156",Zm00001d013156))+geom_beeswarm(aes("Zm00001d006775",Zm00001d006775))+labs(title = "geom_beeswarm")+background
data<-melt(fpkm2[c(1:10),c(1:10,401)],id.vars = "tissue")
names(data)<-c("tissue","gene","FPKM")
tile_plot<-ggplot(data)+geom_tile(aes(gene,tissue,fill=log2(FPKM+1)))+labs(title = "geom_tile")+scale_fill_gradient2(high="red",low="blue",mid="white",midpoint = 6)+background
grid.arrange(density_plot,hist_plot,point_plot,bar_plot,box_plot,violin_plot,bee_plot,tile_plot,nrow=2)






###4 .第一部分画图
fpkm<-read.table("data_for_JCB.txt",header=T,sep = "\t",row.names = "gene")
# 将fpkm 里面大于5000，修改成5000  ;treeheight_row 修改树高度
# My_code
fpkm[fpkm>5000]<-5000
pheatmap(log2(fpkm+1),color = colorRampPalette(c("blue", "white", "red"))(500),
         show_colnames = F, show_rownames = F,clustering_method = "ward.D",
         treeheight_row = 50,
         treeheight_col = 50,border_color = FALSE)


################第一部分基础图形#############
pheatmap(log2(fpkm+1),color = colorRampPalette(c("blue", "white", "red"))(500),
         show_colnames = F, show_rownames = F,clustering_method = "ward.D",
         legend = T, cluster_cols=T,cluster_rows=T,treeheight_row = 50,
         treeheight_col = 50,border_color = FALSE)


head(fpkm[,c(1:10)])
fpkm2<-as.data.frame(t(fpkm))
pheatmap(log2(fpkm2+1),color = colorRampPalette(c("blue", "white", "red"))(500),
         show_colnames = F, show_rownames = F,clustering_method = "ward.D",
         legend = T, cluster_cols=T,cluster_rows=T,treeheight_row = 50,
         treeheight_col = 50,border_color = FALSE)
head(fpkm2[,c(1:3)])



### 5.密度图
library(ggsci)   
# ggsci https://www.jianshu.com/p/71fc7e2561c4
## 一维密度图
head(fpkm2[,c(1:3)])
ggplot(fpkm2)+
  geom_density(aes(Zm00001d013156),fill="grey60",alpha=0.5)+
  labs(title = "geom_density")+
  theme_bw()

# fill 放到aes() 里面。
ggplot(fpkm2)+
  geom_density(aes(Zm00001d013156,fill="red"),alpha=0.5)+
  geom_density(aes(Zm00001d006775,fill="blue"),alpha=0.5)+
  labs(title = "geom_density")+theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+scale_fill_nejm()
  

  

## 二维密度图

ggplot(fpkm2)+
  geom_density2d(aes(Zm00001d013156,Zm00001d011180))+
  labs(title = "geom_density2d")+
  theme_cowplot()

## 热度
ggplot(fpkm2,aes(Zm00001d013156,Zm00001d011180))+
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  geom_point()+geom_smooth(method = "lm")+
  scale_fill_gradientn(colours=c("white","SkyBlue","yellow","red","DarkRed"))+
  labs(title = "stat_density2d")+
  theme_bw()


### 6.柱形图
ggplot(fpkm2)+
  geom_histogram(aes(Zm00001d006775),bins = 30)+
  labs(title = "geom_histogram")+
  theme_bw()


### 7.点图

fpkm2$tissue<-row.names(fpkm2)
head(fpkm2[,c(1:3)])
ggplot(fpkm2)+
  geom_point(aes(tissue,Zm00001d006775))+
  labs(title = "geom_point")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


head(fpkm2[c(30:70),c(1:3)])

ggplot(fpkm2[c(30:70),])+
  geom_point(aes(tissue,Zm00001d006775))+
  labs(title = "geom_point")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))


###8.条形图

ggplot(fpkm2)+
  geom_bar(aes(tissue,Zm00001d006775),stat="identity")+
  labs(title = "geom_bar")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())


###9.条形图极坐标

pie<-data.frame(x=c("a","b","c","d","e"),y=c(12,35,28,9,16))
ggplot(pie)+
  geom_bar(aes(x,y,fill=x),stat="identity")+
  scale_fill_manual(values = colorRampPalette(c("red","purple","blue"))(5))+
  # scale_fill_manual(values = c("red","orange","green","blue","purple"))
  labs(title = "geom_bar")+
  coord_polar(theta = "y") #coord_polar(theta = "x") 
  +theme_cowplot() 
  
# 饼图
ggplot(pie)+
  geom_bar(aes("",y,fill=x),stat="identity")+
  scale_fill_manual(values = c("red","orange","green","blue","purple"))+
  labs(title = "geom_bar")+
  theme_bw() +#coord_polar(theta = "x")
  coord_polar(theta = "y")



### 8.误差线

df <- data.frame(id = c("a", "b", "c", "d"), mid = c(1,5,3,4),
                 upper = c(1.1, 5.3, 3.3, 4.2),lower = c(0.8, 4.6, 2.4, 3.6))

ggplot(df)+
  geom_errorbar(aes(x=id,ymin=lower,ymax=upper),width=0.5,size=1)+
  geom_bar(aes(x=id,y=mid),stat="identity",width=0.8)+
  labs(title = "geom_errorbar")+
  theme_bw()


###9.箱线图

ggplot(fpkm2)+
  geom_boxplot(aes("Zm00001d013156",Zm00001d013156),fill="grey50")+
  labs(title = "geom_boxplot")+
  ylab("")+xlab("")+
  theme_bw()

data<-fpkm2[,c("Zm00001d013156","Zm00001d006775")]
d<-melt(data)

ggplot(d)+
  geom_boxplot(aes(variable,value,fill=variable))+
  scale_fill_manual(values = c("red","blue"))+
  labs(title = "geom_boxplot")+
  ylab("")+xlab("")+
  theme_bw()


### 10.小提琴图
ggplot(fpkm2)+
  geom_violin(aes("Zm00001d013156",Zm00001d013156),fill="grey50")+
  labs(title = "geom_violin")+
  ylab("")+xlab("")+
  theme_bw()

ggplot(d)+
  geom_violin(aes(variable,value,fill=variable))+
  scale_fill_manual(values = c("red","blue"))+
  labs(title = "geom_violin")+
  ylab("")+xlab("")+
  theme_bw()+theme(legend.position = "None")

### 11.蜜蜂图
# stackdir:点图居中对齐
# 使用binaxis='y'选项，则数据点沿y轴进行堆叠，并沿着x轴分组 https://www.jianshu.com/p/b54f7013e67f
ggplot(fpkm2)+
  geom_dotplot(aes("Zm00001d013156",Zm00001d013156),binaxis='y',stackdir='centerwhole',dotsize=0.2)+
  geom_beeswarm(aes("Zm00001d013156",Zm00001d013156),color="red")+
  labs(title = "geom_beeswarm")+
  ylab("")+xlab("")+
  theme_bw()

ggplot(d)+
  geom_beeswarm(aes(variable,value,color=variable))+
  scale_color_manual(values = c("red","blue"))+
  labs(title = "geom_beeswarm")+
  ylab("")+xlab("")+
  theme_bw()


### 12 热图

library(reshape2)
data<-melt(fpkm2[c(1:10),c(1:10,401)],id.vars = "tissue")
names(data)<-c("tissue","gene","FPKM")
ggplot(data)+
  geom_tile(aes(gene,tissue,fill=log2(FPKM+1)))+
  labs(title = "geom_tile")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.3))+
  scale_fill_gradient2(high="red",low="blue",mid="white",midpoint = 5)





################第二部分分面#############

fpkm<-read.table("data_for_JCB.txt",header=T,sep = "\t",row.names = "gene")
fpkm3<-as.data.frame(t(fpkm))[,c(11,157,100,301:303)]
fpkm3$tissue<-row.names(fpkm3)
head(fpkm3)

fpkm3<-melt(fpkm3,id.vars = "tissue")
names(fpkm3)<-c("tissue","gene","FPKM")
head(fpkm3)

tissue_classify<-read.table("sampleID_tissue.txt",sep = "\t")
names(tissue_classify)<-c("tissue","tissue_class")
fpkm3<-merge(fpkm3,tissue_classify,by="tissue")
head(fpkm3)

gene_classify<-data.frame("gene"=c("Zm00001d016256","Zm00001d042346","Zm00001d004380",
                                   "Zm00001d011180","Zm00001d043387","Zm00001d014656"),
                          "gene_class"=rep(c("SPGs","HKGs"), each = 3))
fpkm3<-merge(fpkm3,gene_classify,by="gene")
head(fpkm3)






################基因在不同组织中的表达情况#####################

ggplot(fpkm3)+
  geom_boxplot(aes(gene,log2(FPKM+1),fill=gene_class))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

# 分面 facet_wrap
ggplot(fpkm3)+
  geom_boxplot(aes(gene,log2(FPKM+1),fill=gene_class))+
  facet_wrap(.~tissue_class,nrow = 1)+
  labs("facet_wrap")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))


ggplot(fpkm3)+
  geom_boxplot(aes(gene,log2(FPKM+1),fill=gene_class))+
  facet_wrap(gene_class~tissue_class,nrow = 2)+
  #+facet_wrap(gene_class~tissue_class,scales = "free",nrow = 2)
  labs("facet_wrap")+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

# 分面 facet_grid

(ggplot(fpkm3)
  +geom_boxplot(aes(gene,log2(FPKM+1),fill=gene_class))
  +facet_grid(.~tissue_class)
  +labs("facet_grid")
  +theme_bw()
  +theme(axis.text.x = element_text(angle = 45,hjust=1))
)

(ggplot(fpkm3)
  +geom_boxplot(aes(gene,log2(FPKM+1),fill=gene_class))
  +facet_grid(gene_class~tissue_class)
  #+facet_grid(gene_class~tissue_class,space = "free",scales = "free")
  +labs("facet_grid")
  +theme_bw()
  +theme(axis.text.x = element_text(angle = 45,hjust=1))
)


################第三部分数据图形组合#############
####
head(fpkm3)

ggplot(fpkm3[fpkm3$gene=="Zm00001d014656",],aes(FPKM,..density..))+
  geom_histogram(bins = 20)+
  geom_line(stat="density",size=1,color="blue")+
  geom_density(fill="red",alpha=0.3,color=NA)+
  geom_vline(xintercept = 15,lty=2)+theme_bw()





################geom_boxplot & geom_violin###############

(ggplot(fpkm3,aes(gene,log2(FPKM+1)))
 +geom_violin(aes(fill=gene_class))
 +geom_beeswarm(size=1)
 +geom_boxplot(width=0.1,outlier.colour = NA)
 +geom_hline(yintercept = 2,lty=2)
 +theme_bw()
 +theme(axis.text.x = element_text(angle = 45,hjust=1))
)



#################grid.arrange############################
grid.arrange(P1,P2,P3,P4,P5,P6,
             nrow=2,heights=c(1,1.5),widths=c(1,1,1.5))


grid.arrange(ggplot(),P1,P2,P3,
             nrow=2,heights=c(1,1.5),widths=c(1,1.5))


grid.arrange(ggplot(),P2,ggplot(),
             P3+coord_flip()+scale_y_reverse(),P8,P5+coord_flip(),
             ggplot(),P4,ggplot(),
             nrow=3,heights=c(1,2,1),widths=c(1,2,1))






################第四部分样式精修#############
(ggplot(fpkm3[fpkm3$gene=="Zm00001d014656",],aes(FPKM,..density..))
 +geom_histogram(bins = 30,color="black",fill="grey50")
 +geom_line(stat="density",size=1,color="blue")
 +xlab("FPKMs of Zm00001d014656")+ylab("Density")
 +theme_bw()
 +theme(panel.border = element_rect(colour = "black",size=2))
 +theme(panel.grid = element_line(color = "purple",linetype = 2))## 对所有网格线调整
 +theme(panel.grid.major = element_line(colour = "cyan"))        ## 对大网格线颜色调整
 +theme(panel.grid.minor.x = element_blank())                    ## 对小网格线的纵线进行删除
 +theme(axis.line = element_line(size=2,linetype = 2))           ## 坐标轴线条调整
 +theme(axis.title = element_text(size = 14,colour = "red"))     ## 对坐标轴title进行调整，大小颜色
 +theme(axis.text = element_text(size = 12,colour = "blue"))     ## 对刻度上字进行调整
 +theme(axis.text.x = element_text(angle = 90,vjust = 0.5))      ## 对刻度上x 轴字进行调整
 +theme(axis.ticks = element_line(size=2,color="orange",linetype = 3,arrow = arrow()))   ## arrow 刻度点改成弧形
)






(ggplot(fpkm3[which(fpkm3$FPKM>0),],aes(gene,FPKM))
  +geom_violin(aes(fill=gene_class,color=gene_class))
  +geom_beeswarm(size=1)
  +scale_fill_manual(values = c("grey40","grey70"),labels=c("HouseKeeping","TissueSpecific"))
  +scale_color_manual(values = c("red","blue"))
  +scale_y_log10(breaks=c(0.1,10,100),labels=c("A","D","M"))
  +scale_x_discrete(labels = c("G1","G2","G3","G4","G5","G6"))
  +xlab("Genes")+ylab("FPKM")
  #+geom_boxplot(width=0.1,outlier.colour = NA)
  +labs(fill="Gene Type",color="HKGs & SPGs")
  +theme_bw()
  +theme(axis.text = element_text(size = 12,colour = "black"))
  +theme(axis.title = element_text(size = 14,colour = "black"))
 # +theme(legend.position = "right")
  +theme(legend.justification = "bottom")                         ## 调整图例到右下角
  +theme(legend.position = c(0.8,0.3))                            ##可以直接定位坐标
  +theme(legend.title = element_text(colour = "black",face = "bold"))   ## 图例title 加粗，黑色
  
  
  
  
###################background########################
background<-(
    theme_bw()
    +theme(panel.grid   = element_blank())
    +theme(panel.border = element_rect(size = 0.5))
    +theme(strip.background = element_blank())
    +theme(strip.text   = element_text(size = 12, color = "black"))
    +theme(plot.title   = element_text(size = 15, color = "black", hjust = 0.5, vjust=0.5))
    +theme(axis.ticks   = element_line(size = 0.5,color = "black"))
    +theme(axis.text    = element_text(size = 12, color = "black", angle = 00))
    +theme(axis.title   = element_text(size = 14, color = "black"))
    +theme(legend.text  = element_text(size = 10, color = "black"))
    +theme(legend.title = element_text(size = 10, color = "black" ))
  )

  
)



################第五部分显著性检验#############
library(ggpubr)

head(fpkm3)
compare<-list(c("Zm00001d004380","Zm00001d011180"),
              c("Zm00001d011180","Zm00001d014656"),
              c("Zm00001d016256","Zm00001d014656"))
(ggplot(fpkm3,aes(gene,log2(FPKM+1)))
  +geom_violin(aes(fill=gene_class))
  +geom_beeswarm(size=1)
  +geom_boxplot(width=0.1,outlier.colour = NA)
  +stat_compare_means(comparisons = compare, method = "t.test", label = "p.signif",method.args = list(alternative = "greater"))###wilcox.test
  +theme_bw()
  +theme(axis.text.x = element_text(angle = 30,hjust=1))
)

