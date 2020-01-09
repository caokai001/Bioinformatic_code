###################
## time  : 2020年1月9日20:57:19
## author: kcao
## aim   : overlap riskSNP.bed with ICGC somatic snp
###################

# check String input (sky blue)
if [[ -z $1  ||  $1 == "-h" ]];then
	echo -e "\033[36m Usage:\n\tbash 3.riskSNP_overlap_somatic_SNP.sh Query_rsSNP_file search_base(pattern) \033[0m"
	echo -e "\033[36m Query_rsSNP_file :  
	13	19917046	19917046	rs15
	8       134184336	134184336	rs652          
 search_base(pattern) : 
	simple_somatic_mutation.open.PRAD

										\033[0m"
										echo $1
	exit
fi



# Query_rsSNP_file="3.quary_riskSNP.bed"
# search_base="simple_somatic_mutation"

Query_rsSNP_file=$1
search_base=$2


# chrname=13
# start=19917046
# rsSNP_ID=rs15
###############################

cd /public/home/kcao/Desktop/2019-GR-ERG-AR/IV_check_crisper_rs4919742
# 对每一个rsSNP 遍历

cat $Query_rsSNP_file|while read line;do 
array=($line)
# echo ${array[*]} 
chrname=${array[0]}
start=${array[1]}
rsSNP_ID=${array[3]}

	# 对特定snp 检索不同的ICGC数据库
	ls |grep $search_base |while read id ;
	do cat simple_somatic_mutation.open.PRAD-CA.tsv |
	awk -v chrname=$chrname -v start=$start -v id=$id -v rsSNP_ID=$rsSNP_ID \
	'BEGIN{OFS=FS="\t"}
	($9==chrname && $10==start){print $9,$10,rsSNP_ID,id;exit}'  ;

	done
	
done

