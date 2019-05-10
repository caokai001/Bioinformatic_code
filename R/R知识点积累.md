# R 知识点积累
>1.读取文件 [readr--library(tidyverse)](https://blog.csdn.net/weixin_38423453/article/details/82956517)
 - 分隔符读入：read_csv(), read_csv2(), read_tsv(), read_delim()
 - 空格分隔读入：read_fwf(), read_table()
 - log文件读入：read_log()

>2.read.table()  中出现列名为`X273.eaf`
- 因为R已自动转换了名称，设置为
-df <- read.csv("mydata.csv", check.names=FALSE)
