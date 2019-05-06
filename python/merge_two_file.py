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
