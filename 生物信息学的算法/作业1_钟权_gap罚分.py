import numpy as np
import time

PATH = 'C:\\Users\\zhong\\Documents\\Tencent Files\\1148980644\\FileRecv\\ass1'

with open(PATH + '/input3.fa', 'r') as f:
    f.readline()
    seq1 = f.readline().rstrip()
    f.readline()
    seq2 = f.readline().rstrip()
    
eq = 2
uneq = -1 
score_table = {
        'a': {'a':  eq, 'g': uneq, 'c': uneq, 't': uneq},
        'g': {'a': uneq, 'g':  eq, 'c': uneq, 't': uneq},
        'c': {'a': uneq, 'g': uneq, 'c':  eq, 't': uneq},
        't': {'a': uneq, 'g': uneq, 'c': uneq, 't':  eq}
        }


def fetch_path(seq1, seq2, d=0, s=-1):
    deduce = d + s
    inf = -1000
    len_col = len(seq2) + 1
    len_row = len(seq1) + 1
    m = [0] * len_col
    x = [inf] * len_col
    y = [inf] * len_col
    transition = [0] * (len(seq1) * len(seq2))
    for col in range(1, len_col):
        m[col] = d + col * s
    
    index = 0
    for i in range(1, len_row):
        last_m = m[0]
        m[0] = d + i * s
        for j in range(1, len_col):
            cur_m = m[j]
            x[j] = max((x[j-1] + s, m[j-1] + deduce))
            y[j] = max((y[j] + s, m[j] + deduce))
            
            score = score_table[seq1[i - 1]][seq2[j - 1]]
            max_i, max_value = 0, - 100
            for _i, value in enumerate((last_m + score, x[j], y[j])):
                if value > max_value:
                    max_i = _i
                    max_value = value

            m[j] = max_value
            transition[index] = max_i
            index += 1
            last_m = cur_m
            
    return max_value, transition




def global_align(seq1, seq2):
    st = time.time()
    score, transition = fetch_path(seq1, seq2)
    print(time.time() - st)
    i = len(seq1) - 1
    j = len(seq2) - 1
    align_seq1 = []
    align_seq2 = [] 
    gap = '-'
    while not (i < 0 or j < 0):
        state = transition[i * len(seq2) + j]
        if state == 0:
            align_seq1.append(seq1[i])
            align_seq2.append(seq2[j])
            i -= 1
            j -= 1
        elif state == 1:
            align_seq1.append(gap)
            align_seq2.append(seq2[j])
            j -= 1
        else:
            align_seq1.append(seq1[i])
            align_seq2.append(gap)
            i -= 1
    
    if i < 0:
        align_seq2.extend(list(seq2[:j+1]))
        align_seq1.extend([gap]*(j+1))
    elif j < 0:
        align_seq1.extend(list(seq1[:i+1]))
        align_seq2.extend([gap]*(i+1))

    return ''.join(align_seq1[::-1]), ''.join(align_seq2[::-1])
