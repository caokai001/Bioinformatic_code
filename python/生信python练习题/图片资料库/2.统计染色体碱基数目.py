# 多行注释与取消多行注释

# 多行注释：选定内容，先CTRL+K，然后CTRL+C 
# 取消注释：选定内容，先CTRL+K，然后CTRL+U


import sys
import collections
#args=sys.argv
#filename=args[1]
count_ATCG=collections.OrderedDict() #构建有顺序的字典
filename="./二题.txt"
Bases=["a","t","g","c","n"]
for line in open(filename):
    line=line.strip()
    if line.startswith(">"):
        chr_id=line.strip(">").strip("\t")
        count_ATCG[chr_id]={}
        continue
    line=line.lower()
    for word in line:
        if word not in count_ATCG[chr_id]:
            count_ATCG[chr_id][word]=1
        count_ATCG[chr_id][word]+=1
    print(line)
print(count_ATCG)



def print_information(**count_ATCG):
    for chr in count_ATCG.keys():
        print(chr)
        for base,number in count_ATCG[chr].items():
            print(base,":",number)
        S=sum([x for x in count_ATCG[chr].values()])  ###碱基总数
        print("sum_base:",S)
        if "n" not in count_ATCG[chr]:
            print("n%:0.00%")
        else:
            print("n%:{:.2f}%".format(count_ATCG[chr]["n"]/S*100))
        print("*"*10)

print_information(**count_ATCG)
# >>> print_information(**count_ATCG)
# chr1
# a : 9
# g : 7
# t : 6
# c : 3
# n : 5
# sum_base: 30
# n%:16.67%
# **********
# chr2
# t : 5
# a : 4
# c : 5
# g : 5
# sum_base: 19
# n%:0.00%
# **********
