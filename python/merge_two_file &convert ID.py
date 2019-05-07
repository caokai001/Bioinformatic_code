##差异peak与count里面的区间对应
##运行目录：/public/home/xdyu/kcao/ChIP-seq-BWA-over/8.DESeq2/05.results
import glob
import pandas as pd
deseq2_PATH="/public/home/xdyu/kcao/ChIP-seq-BWA-over/8.DESeq2/04.deseq2"
counts=glob.glob(r"*.count.txt")  ##glob 搜索匹配的路径
for i in range(len(counts)):
	tmp=counts[i][:-10] #'41-EAF1'
	file_c=pd.read_csv(tmp+".count.txt",header=0,sep="\t")
	file_c=file_c[["Geneid","Chr","Start","End","Strand","Length"]]  ###选取某些列合并
	file_d=pd.read_csv(deseq2_PATH+"/"+tmp+".deseq2.output.csv",header=0,sep=",")
	R=pd.merge(file_d,file_c,how="left",left_on="Row.names",right_on="Geneid")
	R.to_csv("/public/home/xdyu/kcao/ChIP-seq-BWA-over/8.DESeq2/05.results/"+tmp+'.result.csv',sep=",")

###Convert gene ID
##python3 接受参数：usage python 脚本名 输入文件 
##refrence https://www.biostars.org/p/22/
##在线工具需要联网
import glob
import pandas as pd
import sys
import mygene
mg = mygene.MyGeneInfo()

A=[]
with open(sys.argv[1],"r") as f:
    for line in f.readlines():
        A.append(line.strip())
##A=A[1:5]
mg = mygene.MyGeneInfo()
B=mg.querymany(A,scopes="ensembltranscript",fields=[ "entrezgene", "symbol","name","alias"],species="human",as_dataframe=True)
B.to_csv(str(sys.argv[1])+".convert.csv")

运行：
for i in *.genelist;do (nohup python gene_convert.py $i &);done



###爬取mygene  https://docs.mygene.info/en/latest/doc/query_service.html#scopes?tdsourcetag=s_pcqq_aiomsg

In [1255]: def fetch(genes): 
      ...:     from urllib.request import urlopen 
      ...:     from urllib.parse import urlencode 
      ...:     import json 
      ...:     data=",".join(genes) 
      ...:     data=urlencode({"q":data,'scope':'ensembl.transcript'}) 
      ...:     req=urlopen("http://mygene.info/v3/query",data=data.encode()) 
      ...:     return json.loads(req.read().decode()) 
In [1255]: genes=['ENSMUST174625'] 
In [1255]: dic=fetch(genes)
