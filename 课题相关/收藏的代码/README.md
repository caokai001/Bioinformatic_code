> 备份课题代码 

1. bedtools_run_AR_total.sh : 包含bedtools 统计及其画图的全部脚本<br>
2. Convert_snp2bed_hg19.R ： 将GWAS snp （无 header）转换成hg19 坐标
3. `3.riskSNP_overlap_somatic_SNP.sh` : 提取查看riskSNP 是否存在于ICGC数据中.
```
 Usage:
	bash 3.riskSNP_overlap_somatic_SNP.sh Query_rsSNP_file search_base(pattern) 
 Query_rsSNP_file :  
	13	19917046	19917046	rs15
	8       134184336	134184336	rs652          
 search_base(pattern) : 
	simple_somatic_mutation.open.PRAD
```
4. merge_mutation_Enrichment_Score_plot.R ： 获得了每一个TF后，对突变得分画图<br>
5. test_delected.R : 两个转录因子交集区域，突变是否更加显著.<br>
6. model_fit_by_R.R : 使用`fitdistrplus`进行模型拟合
