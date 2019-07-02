## 对p值校正

import numpy as np
import pandas as pd
##使用pd.Dataframe 更加方便

def p_adjust_bh(p):
    """BH"""
    a = ["G"] * len(p)
    b = np.arange(1, len(p)+1)
    c = [str(i) + str(j) for i, j in zip(a, b)]
    P = pd.Series(p, index=c)

    stat = np.arange(1, len(P) + 1) * 0.05 / len(P)
    return [P.values.tolist().index(i)+1 for i, j in zip(P[P.argsort().values], stat) if i<j]

p_adjust_bh(p)

################
#p
#Out[131]: [0.053, 0.001, 0.045, 0.03, 0.02, 0.01]
#p_adjust_bh(p)
#Out[132]: [2, 6, 5, 4]
