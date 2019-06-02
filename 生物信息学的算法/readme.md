## 
>## 1. logging 日志
```
import logging
logger = logging.getLogger('L01_1_global_alignment')
logging.basicConfig(level=logging.INFO)
logger.info('########## global alignment ##########')
logger.info('str_one:%r',str_one)
logger.info('str_two:%r',str_two)

```
>## 2.timeit 计时
```
start = timeit.default_timer()
end = timeit.default_timer()
logger.info('running time:%r',str(end-start))
```

>## 3.np.pad ：数组填充功能
```
# Matrix 里填充一行和最后一行，第一列与最后两列
import numpy as np
Matrix = np.arange(1,7).reshape(2,3)    #原始输入数组A
np.pad(Matrix,((1,1),(1,2)),'constant',constant_values = (0,0)) 
```
