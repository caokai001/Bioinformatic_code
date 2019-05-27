# 多行注释
# 选中代码后 快捷键 Ctrl + /
#
# 单行注释
# 选中代码或者光标停留在该行，然后使用快捷键 Ctrl + /
#
# 多行代码缩进
# 选中代码后，快捷键Tab
#
# 多行代码取消缩进
# 选中代码后，快捷键shift + Tab


import shutil
import os
import glob
import pandas as pd
####移动文件shutil ;glob 搜索
os.mkdir("test_data")
files=glob.glob("PubMed*")
path=os.getcwd()
for i in files:
    shutil.move(path+"\\"+i,path+"\\"+test_data)
####
import collections
count_gene=collections.OrderedDict()
file_gtf=r"C:\Users\16926\Desktop\CTCF\test\hg38_chr22.gtf"

##########
with open (file_gtf) as fh:
   for line in fh:
      lineL=line.strip().split("\t")
      chr_num="chr"+lineL[0]
      type_info=lineL[2]
      if type_info=="gene":
         All_gene="all"+"_"+type_info
         if chr_num not in count_gene:
            count_gene[chr_num]={}
            count_gene[chr_num][All_gene]=0
         count_gene[chr_num][All_gene]+=1    #计算all gene个数
         descrip_info=lineL[8]
         descrip_info=descrip_info[:-1]      #主要是有的biotype就是最后一部分了，这时候后面得到的b_type包含“;”,所以一开始就把最后的";"去掉
         descrip_infoL=descrip_info.split("; ")
         biotype=descrip_infoL[4]
         biotypeL=biotype.split(" ")
         b_type=biotypeL[1]
         b_type=eval(b_type)
         if b_type not in count_gene[chr_num]:
              count_gene[chr_num][b_type]=0
         count_gene[chr_num][b_type]+=1

           print(count_gene)
###创建新文件 "test.txt",以tab分隔开
with open("test.txt","w") as f:
    for chr_num,countAll in count_gene.items():
       for k,v in countAll.items():
          print(chr_num+"\t"+k+"\t"+str(v),file=f)
###读入文件画图
#pd_data = pd.read_table("test.txt")

'#########################'
R 画图
'#########################'
setwd("C:\\Users\\16926\\Desktop\\研究生\\【研究生】\\课程作业-曹锴\\BIO-TM\\AGAC_training-annotations\\AGAC_training")
library(ggplot2)

mydata=read.table("test.txt")
p <- ggplot(mydata,aes(x=mydata$V2,y=mydata$V3,fill=mydata$V1))+
  geom_bar(position="dodge",stat="identity")  ###stat="identity":x,y 一一对应
p<-p+xlab("feature")+ylab("count")+labs(fill="染色体")
p+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
#########################
参考：
#######################
#stat_ geom_ 和position:https://zhuanlan.zhihu.com/p/29553873
#ggplot2 画bar图参考:https://blog.csdn.net/RobertChenGuangzhi/article/details/48268099
# 坐标轴倾斜45*，hjust/vjust 向左上方移动

