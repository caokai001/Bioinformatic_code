> 主要的画图代码
```
#!/bin/R
library("tidyverse")
library(ggrepel)

# setwd("/public/home/kcao/Desktop/2018_NC_GS")

# 画图
data=read_delim("merge_mutation_Enrichment_Score",col_names=FALSE,delim=" ")
colnames(data)[1:3]=c("TF_name","sep","score")
data <-data %>% mutate(color=case_when(score >3 ~ "high",2<score & score <3 ~"middle",TRUE ~"low" ))
# str_replace()
data <-data %>%mutate(label=ifelse(score >3,str_replace(data$TF_name,"_MA.*.meme",""),NA))
data$color=factor(data$color,levels=c("high","middle","low"))
P_1 <-ggplot(data,aes(x=TF_name,score))+
	geom_point(aes(col=data$color))+
	# theme(axis.text.x = element_text(angle = 90, hjust = 1))  # 旋转label 方向
	theme_classic()+
	ylim(0.5,6)+
	theme(#axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
	scale_colour_manual(values = c("red","orange", "black"),
						name="Mutation_score",
						labels=c("score >3", "2<score<3","1<score<2" ))+
	# geom_text(label=data$label,check_overlap = TRUE,size = 2,hjust = 0, nudge_x = 0.05)+
	#geom_label(label=data$label)+
	geom_label_repel(label=data$label)+
	ggtitle("2018 GS & 683 Jasper motif")
	
ggsave("merge_mutation_Enrichment_Score.pdf",P_1)


```
![2019年12月18日19:01:39](https://upload-images.jianshu.io/upload_images/9589088-9cefe3b1d9b9fa55.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

