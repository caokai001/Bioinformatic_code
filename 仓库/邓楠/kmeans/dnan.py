
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pylab import * 
mpl.rcParams['font.sans-serif'] = ['SimHei']  ##因为matplotlib库里没有中文字体,所以这个可以解决坐标轴不能显示中文问题
from sklearn import preprocessing
from mpl_toolkits.mplot3d import Axes3D
data = pd.read_csv(r"C:\Users\86182\Desktop\K-means.txt",sep='\t',encoding='gbk')#导入数据
#Z-score标准化
train_x = data[['2019国际排名','2018世界杯','2015亚洲杯']]
scaled_x= preprocessing.scale(train_x)
#k-means做聚类
from sklearn.cluster import KMeans
k= 3
model= KMeans(n_clusters= k, init='k-means++', n_init=10, max_iter= 300)
clf= model.fit(scaled_x)
pre_y= clf.predict(scaled_x)
#concat合并
result= pd.concat([data, pd.DataFrame(pre_y)], axis=1)
result
####方法一画图
for i in range(0,10):
    if result[0][i]==0:
        plt.scatter(result["国家"][i],result[0][i],color='red')
    if result[0][i]==1:
        plt.scatter(result["国家"][i],result[0][i],color='black')
    if result[0][i]==2:
        plt.scatter(result["国家"][i],result[0][i],color='blue')
plt.show()
#####方法二画图
plt.savefig('my_fig.png')
plt.figure(figsize=(8,6))
ax = plt.subplot(projection = '3d')
ax.set_xlabel('2019-International-ranking',color = 'orange',fontsize = 16)
ax.set_ylabel('2018-world-cup',color = 'orange',fontsize = 16)
ax.set_zlabel('2015-asia-cup',color = 'orange',fontsize = 16)

# 绘制3d空间的点
ax.scatter3D(data['2019国际排名'],data['2018世界杯'],data['2015亚洲杯'],c =pre_y,s=90,alpha = 1)