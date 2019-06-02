'''
-*- coding: utf-8 -*-
@Author  : kcao
@Time    : 2019/6/2 10:30
@Software: PyCharm
@File    : L01_1.py
'''
import logging
import numpy as np
import pandas as pd
import timeit  ##timer

start = timeit.default_timer()
logger = logging.getLogger('L01_1_global_alignment')
logging.basicConfig(level=logging.INFO)
##punish_matrix
A=np.ones([5,5],dtype=np.int16)*-1+np.eye(5,dtype=np.int16)*3
A[0,0]=0
name_list=["_","A","C","G","T"]
punish_matrix=pd.DataFrame(A,columns=name_list,index=name_list)  ##创建标签
#    _  A  C  G  T
# _  0 -1 -1 -1 -1
# A -1  2 -1 -1 -1
# C -1 -1  2 -1 -1
# G -1 -1 -1  2 -1
# T -1 -1 -1 -1  2
str_one="ACAATCC"
str_two="AGCATGC"

logger.info('########## global alignment ##########')
logger.info('str_one:%r',str_one)
logger.info('str_two:%r',str_two)
##score_matrix
##初始化矩阵
logger.info('step :initial matrix')
index=["_"]+[i for i in str_one]
columns=["_"]+[i for i in str_two]
score_matrix=pd.DataFrame(np.zeros([len(str_one)+1,len(str_two)+1]),index=index,columns=columns)

##fill the score_matrix
logger.info('step :fill the score_matrix')
punish=-1
for i in range(len(str_one)+1):
    for j in range(len(str_two)+1):
        if i==0 or j ==0:
            score_matrix.iloc[i,j]=0+punish*i+punish*j
        else:
            insert=score_matrix.iloc[i,j-1]+punish
            delect=score_matrix.iloc[i-1,j]+punish
            match=score_matrix.iloc[i-1,j-1]+punish_matrix.loc[str_one[i-1],str_two[j-1]]
            score_matrix.iloc[i, j]=max(insert,delect,match)

##back_tracking the alignment
logger.info('step :back_tracking the alignment')
P=[]
i=len(str_one)
j=len(str_two)
while i or j :
    if (i>0 and j>0 and score_matrix.iloc[i, j]==score_matrix.iloc[i-1,j-1]+punish_matrix.loc[str_one[i-1],str_two[j-1]]):
        P.append([score_matrix.index[i],score_matrix.columns[j]])
        i=i-1
        j=j-1
    elif (j>0 and score_matrix.iloc[i, j]==score_matrix.iloc[i,j-1]+punish):
        P.append(["_",score_matrix.columns[j]])
        j=j-1
    elif (i>0 and score_matrix.iloc[i, j]==score_matrix.iloc[i-1,j]+punish):
        P.append([score_matrix.index[i],"_"])
        i=i-1

##print result
logger.info('step :print result ')
A_str=B_str=""
for i,y in P[::-1]:
    A_str = A_str + i
    B_str = B_str + y
print("str_one : {}\nstr_two : {}".format(A_str,B_str))

##done!
logger.info('done ^-^ ！')

end = timeit.default_timer()
logger.info('running time:%r',str(end-start))






###################################################
# 结果
###################################################

INFO:L01_1_global_alignment:########## global alignment ##########
INFO:L01_1_global_alignment:str_one:'ACAATCC'
INFO:L01_1_global_alignment:str_two:'AGCATGC'
INFO:L01_1_global_alignment:step :initial matrix
INFO:L01_1_global_alignment:step :fill the score_matrix
INFO:L01_1_global_alignment:step :back_tracking the alignment
INFO:L01_1_global_alignment:step :print result 
str_one : A_CATCC
str_two : AGCATGC
INFO:L01_1_global_alignment:done ^-^ ！
INFO:L01_1_global_alignment:running time:'0.031113404366736708'
