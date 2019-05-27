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
with open(file_gtf) as fh:
    for line in fh.readlines():
        line_L=line.strip().split("\t")
        chr_num="chr"+line_L[0]
        type_info=line_L[2]
        if type_info=="gene":
            All_gene="all"+"_"+type_info
            if chr_num not in count_gene:
                count_gene[chr_num] = {}
                count_gene[chr_num][All_gene]=0
            count_gene[chr_num][All_gene]+=1
            descript_info=line_L[8]
            descript_info_L=descript_info.split("; ")
            biotype=descript_info_L[4]
            biotype_L=biotype.split(" ")
            b_type=biotype_L[1]
            b_type=eval(b_type)
            if b_type not in count_gene[chr_num]:
                count_gene[chr_num][b_type]=0
            count_gene[chr_num][b_type]+=1
for chr_num,count_ALL in count_gene.items():
    for k,v in count_ALL.items():
        print(chr_num,k,v)
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
for chr_num,countAll in count_gene.items():
   for k,v in countAll.items():
      print(chr_num,k,v)
