# 参考文章：https://www.jianshu.com/p/e88190b1e1ed
# 小基因组读入字典，问题不大。

# import sys
# import os
#
# args = sys.argv
# ref_file = args[1]
# chr_info = args[2]
# base_info = int(args[3]) - 1
# print(base_info)
# import collections
#
# os.chdir(r"C:\Users\16926\Desktop\项目\课题进展\bam-readcount")
# aDict = collections.OrderedDict()
# # with open("ref_file") as f:
# with open("HHH_ref.txt") as f:
#     for line in f.readlines():
#         line = line.strip()
#         if len(line) == 0:
#             raise IOError("blank row")
#         if line.startswith(">"):
#             tmp = line.lstrip(">")
#             aDict[tmp] = ""
#         else:
#             aDict[tmp] += line
# print(aDict[chr_info], base_info)
# print(aDict[chr_info][base_info])

###
import sys
import os

args = sys.argv
ref_file = args[1]
chr_info = "chr2"
base_up = int(args[3]) - 3
base_down = int(args[3])-1
import collections

os.chdir(r"C:\Users\16926\Desktop\项目\课题进展\bam-readcount")
aDict = collections.OrderedDict()
state = False
n = 0

with open("HHH_ref.txt") as f:
    for line in f.readlines():
        lineL = line.strip()
        if lineL.startswith(">") and lineL[1:] == chr_info:
            state = True

        if state:
            n += len(line)
            if n <= base_up:
                continue
            else:
                print("base_up", lineL[-(n - base_up)])

            if n >= base_down:
                print("base_down", lineL[-(n - base_down)])
            else:
                continue

            break
            
###
运行：$ python chr_base.py HMM_ref.txt chr2 15
base_up A
base_down G
