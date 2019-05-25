#C:\Users\hp\Desktop\学习\python\hzau.python\第一题.py
#参考链接：https://www.jianshu.com/p/8cddd0774b08
#下载链接：wget ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
#tips：考虑了有些基因公用exon,加入字典，进行过滤筛选

import re
import os
filename="CCDS.current.txt"
exon_length=0
aDict={}
i=1
os.chdir(r"c:\Users\hp\Desktop")
with open(filename) as fh:
    for line in fh.readlines():
        line=line.strip()
        if line.startswith("#"):
            continue
        line=line.split("\t")
        exon=line[-2]
        if exon=="-":
            continue
        exonT=re.sub("\[|\]","",exon).split(",")  ##字符串replace
        for exon in exonT:
            exon_name=line[1]+":"+exon
            if exon_name not in aDict:
                aDict[exon_name]=1
                start=int(exon.split("-")[0].strip())
                end=int(exon.split("-")[1].strip())
                exon_length += end-start
print(exon_length)
#line="1	NC_000001.11	VWA1	64856	CCDS27.1	Public	+	1435748	1439786	[1435748-1435820, 1436926-1437483, 1439080-1439786]	Identical"
