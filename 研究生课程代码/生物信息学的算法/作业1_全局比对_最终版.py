# Optimize spatial time and backtracking process to make it suitable for long sequence global alignment
# filling score matrix row by row
# use mid point algorithm to divide the sequence to two substring
# then use Recursive to output the alignment information

# usage :python script_name.py parameter.txt input.fa output.txt

import logging
import os
import sys
import time
import numpy as np
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("parameter_txt",help="parameter file ,while contain score for initiating a gap;\
                                      score for each base insert/delete;alphabet;similarity matrix  ")

parser.add_argument("input_fa",help="input fasta format")
parser.add_argument("output_txt",help="the name of output file")
parser.add_argument('--version', action='version', version='version 1.0')


args=parser.parse_args()
# print(args["input.fa"])
parameter=args.parameter_txt
input_fa=args.input_fa
output_file=args.output_txt

########################################
# Parameter information
logger = logging.getLogger('ss1_global_alignment')
logging.basicConfig(level=logging.INFO)
#os.chdir(r"C:\Users\16926\Desktop\研究生\【研究生】\研究生课程\算法课程\2019-生信算法课\ass1")
#parameter, input_fa, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
global_max_score = None
########################################

# 读取parameter文件
count = 0
similarity_matrix = []
with open(parameter) as f:
    for line in f.readlines():
        line = line.strip().strip("\n")
        if len(line):
            count += 1
            if count == 1:
                init_gap = line.split(";")[0]
            elif count == 2:
                indel = int(line.split(";")[0])
            elif count == 4:
                alphabet = line.split(" ")
            elif count > 5:
                similarity_matrix.append(line.split())

# loading input file
seq = {}  ##save seq
with open(input_fa) as f:
    for line in f.readlines():
        if line.startswith(">"):
            name = line.lstrip(">").strip("\n")
            seq[name] = ""
        else:
            seq[name] += line.strip("\n").strip()

str_one = seq["seq1"]
str_two = seq["seq2"]


# fill and memory the position by row
# punish_matrix score_matrix
punish_score={}
#alphabet.append("-")
#['a', 'c', 'g', 't']
for i in range(len(alphabet)):
    punish_score[alphabet[i]]={}
    for j in range(len(alphabet)):
        punish_score[alphabet[i]][alphabet[j]]=int(similarity_matrix[i][j])
    punish_score[alphabet[i]]["-"] = -1
punish_score["-"]={}
for i in alphabet:
    punish_score["-"][i] = -1

# punish_score = {'-': {'a':-1, 'g': -1, 'c': -1, 't': -1, '-': -1},
#                 'a': {'a': 2, 'g': -1, 'c': -1, 't': -1, '-': -1},
#                 'c': {'a': -1, 'g': -1, 'c': 2, 't': -1, '-': -1},
#                 'g': {'a': -1, 'g': 2, 'c': -1, 't': -1, '-': -1},
#                 't': {'a': -1, 'g': -1, 'c': -1, 't': 2, '-': -1}}


def find_score(string_one, string_two):
    '''
    find_score() function show the max scale of two string:
    sush as:

    str_one="acaatcc"
    str_two="agcatgc"

    find_score(str_one,str_two)=7
    it means The maximum alignment score of str_one and str_two is 7.

    find_score(str_one[0:4],str_two[0:7])=V[len(str_two[0:7]=2

    '''
    n = len(string_one)
    m = len(string_two)
    V = [0] * (m + 1)
    for j in range(1, m + 1):
        V[j] = V[j - 1] + indel
    for i in range(1, n + 1):
        V_prev = V[0]
        V[0] = V_prev + indel
        for j in range(1, m + 1):
            V_cur = V[j]
            V[j] = max(V_prev + punish_score[string_one[i - 1]][string_two[j - 1]],
                       V_cur + indel,
                       V[j - 1] + indel)
            V_prev = V_cur
    return V
# find_score(str_one[0:4],str_two[0:7])

# mid-point_agrithom
def mid_point(string_one, string_two):
    '''

    :param string_one:
    :param string_two:
    :return: the breakpoint position.
    For example:
    (4,4) means that the first four characters of the two strings are mid_point positions.
    '''
    n = len(string_one)
    m = len(string_two)
    i_mid = n // 2

    score1 = find_score(string_one[0:i_mid], string_two)
    score2 = find_score(string_one[i_mid:n][::-1], string_two[::-1])

    max_sum = -100
    max_j = 0
    length = len(score1) - 1
    for j in range(length + 1):
        sums = score1[j] + score2[length - j]
        if sums > max_sum:
            max_j = j
            max_sum = sums

    global global_max_score
    if global_max_score is None:
        global_max_score = max_sum

    return i_mid, max_j
# mid_point("acaatcc","agcatgc")


def back_track(string_one, string_two):
    if not string_two:
        return [string_one, "-" * len(string_one)]

    if (len(string_one) == 1) or (len(string_two) == 1):
        if len(string_one) == 1:
            Direction = 1
            # print("string_two=1"+"#"*20)
            # print(string_one,string_two)
        else:
            string_one, string_two = string_two, string_one
            Direction = -1

        index = np.argmax([punish_score[string_one][ch] for ch in string_two])
        return ['-' * (index) + string_one + '-' * (len(string_two) - index - 1),
                string_two][::Direction]

    else:
        p, q = mid_point(string_one, string_two)
        return back_track(string_one[0:p], string_two[0:q]) + back_track(string_one[p:], string_two[q:])


def global_align(s1, s2):
    array = back_track(s1, s2)
    S1, S2 = (array[i] for i in range(0, len(array), 2)), (array[i] for i in range(1, len(array), 2))
    output_file.writelines("score = {:.1f}\n".format(global_max_score))
    output_file.writelines(">seq1\n" + "".join(S1) + "\n" + ">seq2\n" + "".join(S2))
    #print(">seq1\n" + "".join(S1) + "\n" + ">seq2\n" + "".join(S2))


# main program
if __name__ == '__main__':
    logger.info("     loading argument:\n     %s ,%s ,%s",parameter,input_fa,output_file)
    output_file = open(output_file, "w",encoding='utf8')
    start = time.time()
    global_align(str_one, str_two)
    # print(find_score(str_one, str_two))
    end = time.time()
    logger.info('     running time: %ss',str(end - start))
    logger.info("     score = {:.1f}".format(global_max_score))
    logger.info("     finished !")
    logger.info("    !!please open the output file to see more information!! ")
    output_file.close()
