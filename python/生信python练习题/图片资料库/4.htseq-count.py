import glob
import collections
mydict=collections.OrderedDict()
list_dirs=glob.glob("[0-9].txt")
for i in list_dirs:
   for line in open (i):
      array=line.strip().split("\t")
      if array[0] not in mydict:
           mydict[array[0]]=[array[1]]    #注意array[1]外面的中括号
      else:
           mydict[array[0]].append(array[1])    #字典中，对一个key增加多个value的方法

for i in mydict:
     print("{0}\t{1}".format(i, "\t".join(mydict[i])))
#########
# 结果
#gene1	2	6
#gene2	6	78
#gene3	6	9
#gene4	8	6
#gene7	3	36
