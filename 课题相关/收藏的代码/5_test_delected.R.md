>考虑转录因子之间组合，TF-TF 突变得分分布4个图，参考皮肤癌CTCF_RAD21

输入数据url：https://github.com/caokai001/Bioinformatic_code/blob/master/%E8%AF%BE%E9%A2%98%E7%9B%B8%E5%85%B3/%E5%A4%87%E4%BB%BD%E6%95%B0%E6%8D%AE%E5%A4%B9/test_delected.txt

```
mutation_paired_TF <-file.choose()  # "C:\\Users\\16926\\Desktop\\test_delected.txt"
data <- read.table(mutation_paired_TF,header = FALSE)
colnames(data) <-c("tf_name","co_tf_name","tf_without_score","tf_within_score","co_tf_within_score","co_tf_without_score")

head(data)
data[data$tf_name=="Ahr::Arnt_MA0006.1.meme" & data$co_tf_name=="NRF1_MA0506.1.meme",]
data[data$tf_name=="NRF1_MA0506.1.meme" & data$co_tf_name=="Ahr::Arnt_MA0006.1.meme",]




## 得到候选组合
data_log <- cbind(data[1:2],-log10(data[-1:-2]))
filter_data <- data_log[which(data_log$tf_without_score < data_log$tf_within_score),]
dim(filter_data)

## library(tidyverse)
tmp_data <-tbl_df(filter_data) %>%
  mutate("delta_2_1" =tf_within_score -tf_without_score) %>%
  arrange(desc(delta_2_1)) %>%
  filter(delta_2_1 >2)
View(tmp_data)

## 绘制delta 变化范围
library(scales)
library(ggplot2)
P_1 <-ggplot(tmp_data,aes(delta_2_1)) +
  geom_histogram(bins = 300)+
  # geom_density()+
  scale_x_continuous(trans=log2_trans()) +
  xlab("delta size")+
  scale_y_continuous(breaks = seq(1,10))+
  ggtitle("calculate the difference between TF without and within cofactor")+
  theme(plot.title = element_text(size=12))
ggsave("delta_score_distribution.pdf" ,P_1,width = 24, height = 6)
## 尝试热图绘制
# https://guangchuangyu.github.io/cn/2019/03/pheatmap-ggplotify/
library(pheatmap)
require(ggplotify)

data_heat <- tmp_data[3:6]
rownames(data_heat) <- paste(tmp_data$tf_name,tmp_data$co_tf_name,sep=" || ")
p<-pheatmap(data_heat,scale = "row",show_rownames=F,cluster_cols = F,angle_col=0)
P <-g+ggtitle("  The mutation score of different TF type\n")
ggsave("The_mutation_score_of_different_TF_type.pdf" ,P,width = 6, height = 10)

```

例：两个转录因子交互区域，突变是否更强<br>
![2020年1月10日22:08:26](https://upload-images.jianshu.io/upload_images/9589088-b4de8ac43b546c61.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

例：输入文件展示<br>
![2020年1月10日22:08:40](https://upload-images.jianshu.io/upload_images/9589088-6a1373195d39c249.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

例：输入数据，用pheatmap 热图聚类<br>
![2020年1月10日22:10:08](https://upload-images.jianshu.io/upload_images/9589088-fcc15e9bdbd5d76a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
