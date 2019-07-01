### 1：boundary 附近的ChIP 分布图
```
$ awk '{if($4==0){print $0}}' prec250.tad |sort -Vk 1 -k2n |cut -f 1,2,3 >PREC.tad.bed
[kcao@comput7 2.bw]$ nohup computeMatrix scale-regions -p 10 -R PREC.tad.bed -S SRR513122.deeptools.bw -b 500000 -a 500000 --regionBodyLength 100000 --skipZeros -o PREC.heatmap.gz&
[kcao@comput7 2.bw]$ plotHeatmap -m merged_TSS.gz -out merged.png
[kcao@comput7 2.bw]$ plotProfile --dpi 720 -m PREC.heatmap.gz --startLabel left_boundary --endLabel right_boundary --samplesLabel PREC_CTCF_profile_boundary  -out PREC_CTCF_profile_TAD.pdf --plotFileFormat pdf&
```
![效果如下，不是很明显富集，不知道原因](https://upload-images.jianshu.io/upload_images/9589088-2a678f68f5b99f25.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
