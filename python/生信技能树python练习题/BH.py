import numpy as np
import pandas as pd


def p_adjust_bh(p):
    """BH"""
    a = ["G"] * len(p)
    b = np.arange(1, len(p)+1)
    c = [str(i) + str(j) for i, j in zip(a, b)]
    P = pd.(p, index=c)

    stat = np.arange(1, len(P) + 1) * 0.05 / len(P)
    return [P.values.tolist().index(i)+1 for i, j in zip(P[P.argsort().values], stat) if i<j]

p_adjust_bh(p)
