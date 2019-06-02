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
