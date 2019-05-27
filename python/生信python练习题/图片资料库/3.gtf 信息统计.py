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
