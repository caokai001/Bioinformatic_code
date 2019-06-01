-*- coding: utf-8 -*-
@Author  : kcao
@Time    : 2019/6/1 13:59
@Software: PyCharm
@File    : L3.rearrange.py
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
# import sys
# argvs=sys.argv

logger = logging.getLogger('L3_arrange.py') ##logging 名字
logging.basicConfig(level=logging.INFO)  # basicConfig是logging提供的简单的配置方法，不用basicConfig则需要手动添加handler


#A = [8, 9, 3, 4, 7, 6, 5, 1, 2, 10, 11]
A=[4, 5, 3, 1, 2]
logger.info('step A:%r',A) ##格式化输出

def bp_number(a):
    '''
    return 返回的函数包括：
    num:breakpoint
    up_strip:上升序列

    '''
    a = a.copy()## 保留A序列
    a.insert(0, 0)
    a.append(max(a) + 1)
    # [0, 8, 9, 3, 4, 7, 6, 5, 1, 2, 10, 11, 12]
    num = []
    # num 储存第几个逗号是breakpoint
    for i in range(1, len(a)):
        if abs(a[i] - a[i - 1]) != 1:
            num.append(i)
    down_strip = []
    ##用断点来计算up
    for i in range(len(num) - 1):
        if a[num[i]] - a[num[i + 1] - 1] >= 0:     # Strip of size 1 is assumed to be decreasing.
            down_strip.append([num[i]-1, num[i + 1] - 2])
    return [num, down_strip]

###
def min_index(b):
    '''
    min_index(bp_number(A)[1])
    返回list 里值最小的index
    b 为降序序列区间index
    '''
    min_index=A.index(max(A))
    for i in bp_number(A)[1]:                       # decreasing 中最小的index
        if A[i[1]]<A[min_index]:
            min_index=i[1]
    return min_index


while len(bp_number(A)[0]) > 0:                     # 当断点为0时，循环终止
    if len(bp_number(A)[1]) > 0:                    # 降序的子序列数目大于0
        P_min = min_index(bp_number(A)[1])
        if A[P_min]==1:
            A=list(reversed(A[0:P_min+1]))+A[P_min+1:]
            logger.info('step A:%r',A)              ##格式化输出
            continue
        P_min_1=A.index(A[P_min]-1)
        if P_min_1<P_min :                     ##P_min_1在左边reversal
            A = A[:P_min_1 + 1] + list(reversed(A[P_min_1+1:P_min+1])) + A[P_min + 1:]
            logger.info('step A:%r',A)              ##格式化输出
        else:
            A=A[:P_min+1]+list(reversed(A[P_min+1:P_min_1+1]))+A[P_min_1+1:]
            logger.info('step A:%r',A)              ##格式化输出
    else:
        x,y=bp_number(A)[0][0:2]
        A=A[:x-1]+list(reversed(A[x-1:y-1]))+A[y-1:]
        logger.info('step A:%r',A)                  ##格式化输出








###########################################################
##结果：
INFO:L3_arrange.py:step A:[4, 5, 3, 1, 2]
INFO:L3_arrange.py:step A:[4, 5, 3, 2, 1]
INFO:L3_arrange.py:step A:[1, 2, 3, 5, 4]
INFO:L3_arrange.py:step A:[1, 2, 3, 4, 5]
