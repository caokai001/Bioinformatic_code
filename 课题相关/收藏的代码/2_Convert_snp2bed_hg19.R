if (!requireNamespace("biomaRt", quietly = TRUE))
  install.packages("biomaRt")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
if (!requireNamespace("naturalsort", quietly = TRUE))
  install.packages("naturalsort")

# load R packages
suppressMessages(library(biomaRt))
suppressMessages(library(tidyverse))
#suppressMessages(library(naturalsort))

# Parsing command line arguments
# Set to "T", the first input parameter, indexed by 1
args=commandArgs(T)


if(length(args)==0){
  cat("Usage:
      Rscript Convert_snp2bed_hg19.R riskSNP_file\n")
  stop("please input args\n")
}else{
  riskSNP_path=args[1]
  cat(paste0(riskSNP_path," is working...\n"))
}


# snp_list <- read.table(file = "https://raw.githubusercontent.com/caokai001/Bioinformatic_code/master/%E8%AF%BE%E9%A2%98%E7%9B%B8%E5%85%B3/%E5%A4%87%E4%BB%BD%E6%95%B0%E6%8D%AE%E5%A4%B9/GWAS_riskSNP_review_rsID.bed",header = F)
snp_list <- read.table(file=riskSNP_path,header = FALSE)

snp_ids <- as.character(snp_list$V1)
#snp_ids = c("rs746229459","rs4907792")
# snp_ids <-head(snp_ids)

human_variation = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp",host="http://grch37.ensembl.org") 

snp_attributes = c("refsnp_id", 
                   "ensembl_gene_stable_id", 
                   "chr_name", 
                   "chrom_start", 
                   "chrom_end", 
                   "ensembl_transcript_stable_id");
snp_locations = getBM(attributes = snp_attributes, 
                      filters = "snp_filter", 
                      values = snp_ids, 
                      mart = human_variation);

# remove repeat item
snp_chrpos <- tbl_df(snp_locations) %>% arrange(refsnp_id,nchar(chr_name)) %>%
  distinct(refsnp_id,.keep_all = T) %>%
  dplyr::select(refsnp_id,chr_name,chrom_start,chrom_end,refsnp_id) 



# View(snp_chrpos)

# add chr to column 2
merged_data <-  merge(snp_list,snp_chrpos,by.x="V1",by.y="refsnp_id",all.x=TRUE) %>%
  mutate(chr_name=paste0("chr",chr_name),genome="hg19") %>%
  dplyr::select(chr_name,chrom_start,chrom_end,V1,genome) %>% 
  rename(rsSNP=V1) %>%
  arrange(as.numeric(str_replace(chr_name,"chr","")),chrom_start)

# save file
filename=paste0(str_split(riskSNP_path,"\\.",simplify = T)[1,1],"_riskSNP_hg19.bed")
write_delim(merged_data,path = filename ,delim = "\t")
cat("Done !!!\n")
