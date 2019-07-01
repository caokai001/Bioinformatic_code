##参考：https://www.jianshu.com/p/af01a6f194b1
##对hg38.tss储存为染色体为key,value为字典的有序字典。每个基因区间作为key,基因名作为value.
## 用区间信息进行判断，是不是落在此基因上。

import wget
import re

aDict = {}
url = "http://www.biotrainee.com/jmzeng/tmp/hg38.tss"
wget.download(url, out=re.sub(".*/", "", url))

for line in open("hg38.tss"):
    lineL = line.strip().split("\t")
    gene_ID = lineL[0]
    chr_name = lineL[1]
    start = lineL[2]
    end = lineL[3]
    if chr_name not in aDict:
        aDict[chr_name] = {}
    aDict[chr_name][start, end] = gene_ID

for line in open("file.txt"):
    lineList = line.strip().split()
    Chr_name = lineList[0]
    Start = int(lineList[1])
    End = int(lineList[2])
    for K1, V1 in aDict[Chr_name].items():
        if int(K1[0]) <= Start <= int(K1[1]) or int(K1[0]) <= Start <= int(K1[1]):
            print(Chr_name, Start, End, V1)
            
            
