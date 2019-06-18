###1. 评估问题（向前算法）
import numpy as np

def HMM(Alignment,value):
    state = ['Sunny', 'Cloudy', 'Rainy']  # 此为天气状态序列
    seaweed = ['Dry', 'Dryish', 'Damp', 'Soggy']  # 此为海草状态序列
    pi = np.mat([0.63, 0.17, 0.20])  # 存储起始概率
    transfer = np.mat([[0.500, 0.250, 0.250],  # 存储状态转移矩阵
                       [0.375, 0.250, 0.375],
                       [0.125, 0.675, 0.200]])
    output = np.mat([[0.60, 0.20, 0.15, 0.05],  # 存储输出矩阵，列分别代表Dry, Dryish, Damp, Soggy
                     [0.25, 0.25, 0.25, 0.25],
                     [0.05, 0.10, 0.35, 0.50]])
    alpha=np.multiply(pi.T,output[:,seaweed.index(Alignment[0])])
    for i in range(1,len(Alignment)):
        tmp=np.dot(transfer.T,alpha)
        alpha=np.multiply(tmp,output[:,seaweed.index(Alignment[i])])
    return sum(alpha)
    
###2.解码问题（Viterbi算法）
X=1,6,2,6
import pandas as pd
data=[0.5,0.5]
index=["F","L"]
Pi=pd.DataFrame({"state":data},index=index)
Alpha=pd.DataFrame({"F":[0.95,0.05],"L":[0.05,0.95]},index=["F","L"])
E_matrix=pd.DataFrame({"F":[1/6]*6,"L":[1/10]*5+[1/2]},index=range(1,7))

##第一步
tmp=np.multiply(E_matrix.loc[X[0],:].values,Pi.values.T)
index={}
for j in range(1,len(X)):
    V = []
    index[j]=[]
    for i in ["F","L"]:
        A=np.multiply(tmp,Alpha.loc[:,i].values)*E_matrix.loc[X[j],i]
        tmp_H = np.max(A)
        index[j].append(np.where(np.argmax(A)==0,"F","L").tolist())
        V.append(tmp_H)
    tmp=V
#print(tmp)
#print(index)
cur_new=np.where(np.argmax(tmp)==0,"F","L").tolist()


str=[cur_new]
for i in range(len(X)-1,0,-1):
    cur_new=index[i]["FL".index(cur_new)]
    str.append(cur_new)
print(str[::-1])
