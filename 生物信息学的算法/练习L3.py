'''
-*- coding: utf-8 -*-
@Author  : kcao
@Time    : 2019/6/1 13:59
@Software: PyCharm
@File    : L3.rearrange.py

tips:
    时间复杂度为O(M*N),长序列耗时长
    空间复杂度为O(M*N)

usage: python L01_1.py parameter1.txt input1.fa output.txt
全局比对
'''


# reversal methods

########################################
#参考课程：L03_genome_rearrangement.pdf
# 实现 翻转的算法。效率大约是最优方法的4倍。所以叫

#  4-approximation algorithm (IV)
# Algorithm simpleApprox
# • while b(π) > 0,
# • if there exist a decreasing strip,
# • we reverse π by ρmin [this reversal reduces b(π) by at least 1];
# • else
# • reverse an increasing strip to create a decreasing strip [b(π) does not change]
########################################



import logging
import os
import numpy as np
import pandas as pd
import timeit  ##timer
import sys

os.chdir(r"C:\Users\16926\Desktop\研究生\【研究生】\研究生课程\算法课程\2019-生信算法课\ass1")

#parameter = sys.argv[1]
#input_fa = sys.argv[2]
#output_file = sys.argv[3]

parameter,input_fa,output_file="parameter2.txt","input2.fa","output.txt"

start = timeit.default_timer()
logger = logging.getLogger('L01_1_global_alignment')
logging.basicConfig(level=logging.INFO)

output_file=open(output_file,"w")

## 读取parameter文件
count=0
similarity_matrix = []
with open(parameter) as f:
    for line in f.readlines():
        line = line.strip().strip("\n")
        if len(line):
            count += 1
            if count == 1:
                init_gap=line.split(";")[0]
            elif count ==2:
                indel=line.split(";")[0]
            elif count ==4:
                alphabet=line.split(" ")
            elif count >5:
                similarity_matrix.append(line.split())



## 读取input 文件
seq={}    ##储存seq
with open(input_fa) as f:
    for line in f.readlines():
        if line.startswith(">"):
            name=line.lstrip(">").strip("\n")
            seq[name] = ""
        else:
            seq[name] += line.strip("\n").strip()


##punish_matrix
## np.pad 矩阵填充
A=np.array(similarity_matrix,dtype="int32")
A=np.pad(A,((1,0),(1,0)),'constant', constant_values=str(indel))
A[0,0]=init_gap

name_list=["-"]+alphabet
punish_matrix=pd.DataFrame(A,columns=name_list,index=name_list)  ##创建标签
#     _   a   c   g   t
# _   0  -1  -1  -1  -1
# a  -1   2  -1  -1  -1
# c  -1  -1   2  -1  -1
# g  -1  -1  -1   2  -1
# t  -1  -1  -1  -1   2



#str_one="TGAAGTC"
#str_two="TAAGGC"
str_one=seq["seq1"]
str_two=seq["seq2"]


logger.info('########## global alignment ##########')
logger.info('str_one:%r',str_one)
logger.info('str_two:%r',str_two)
##score_matrix
##初始化矩阵
logger.info('step :initial matrix')
index=["-"]+[i for i in str_one]
columns=["-"]+[i for i in str_two]
score_matrix=pd.DataFrame(np.zeros([len(str_one)+1,len(str_two)+1]),index=index,columns=columns)

##fill the score_matrix
logger.info('step :fill the score_matrix')

for i in range(len(str_one)+1):
    for j in range(len(str_two)+1):
        if i == 0 or j == 0:
            score_matrix.iloc[i,j] = int(init_gap)+int(indel)*i+int(indel)*j
        else:
            insert = score_matrix.iloc[i,j-1]+int(indel)
            delect = score_matrix.iloc[i-1,j]+int(indel)
            match = score_matrix.iloc[i-1,j-1]+int(punish_matrix.loc[str_one[i-1],str_two[j-1]])
            score_matrix.iloc[i, j] = max(insert,delect,match)

##back_tracking the alignment
logger.info('step :back_tracking the alignment')
P=[]
i=len(str_one)
j=len(str_two)
while i or j :
    if (i>0 and j>0 and score_matrix.iloc[i, j]==score_matrix.iloc[i-1,j-1]+int(punish_matrix.loc[str_one[i-1],str_two[j-1]])):
        P.append([score_matrix.index[i],score_matrix.columns[j]])
        i=i-1
        j=j-1
    elif (j>0 and score_matrix.iloc[i, j]==score_matrix.iloc[i,j-1]+int(indel)):
        P.append(["-",score_matrix.columns[j]])
        j=j-1
    elif (i>0 and score_matrix.iloc[i, j]==score_matrix.iloc[i-1,j]+int(indel)):
        P.append([score_matrix.index[i],"-"])
        i=i-1


##print result
score=0
for x,y in P:
    score+=int(punish_matrix.loc[x,y])


logger.info('step :print result ')
A_str=B_str=""
for i,y in P[::-1]:
    A_str = A_str + i
    B_str = B_str + y
#print("score = {}".format(score))
#print("str_one : {}\nstr_two : {}".format(A_str,B_str))
output_file.write("score = {:.1f}\n>seq1\n{}\n>seq2\n{}".format(score,A_str,B_str))
output_file.close()

print("score = {:.1f}\n>seq1\n{}\n>seq2\n{}".format(score,A_str,B_str))

##done!
logger.info('done ^-^ ！')

end = timeit.default_timer()
logger.info('running time:%r',str(end-start))





#############################################
#结果
############################################
INFO:L01_1_global_alignment:########## global alignment ##########
INFO:L01_1_global_alignment:str_one:'HEAGAWGHEE'
INFO:L01_1_global_alignment:str_two:'PAWHEAE'
INFO:L01_1_global_alignment:step :initial matrix
INFO:L01_1_global_alignment:step :fill the score_matrix
INFO:L01_1_global_alignment:step :back_tracking the alignment
INFO:L01_1_global_alignment:step :print result 
score = -8.0
>seq1
HEAGAWGHEE
>seq2
--P-AWHEAE
INFO:L01_1_global_alignment:done ^-^ ！
INFO:L01_1_global_alignment:running time:'0.04599526172114565'


