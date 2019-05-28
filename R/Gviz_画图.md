### 1.画bam 文件
```
library(Rsamtools)
setwd("C:\\Users\\16926\\Desktop\\CTCF")
library(rtracklayer)
library(Gviz)
indexBam("SRR513122.chr.uniq.bam")
bamFile <- system.file("SRR513122.chr.uniq.bam", package = "Gviz")

Bam_file <- DataTrack(range = "SRR513122.chr.uniq.bam", genome = "hg38",
                     type = "histogram", name = "chr17-CTCF", window = -1, chromosome = "chr17",
                     windowSize = 100,fill = "blue",fill.histogram = "blue",
                     ylim = c(0, 15))
idxTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr17")

axisTrack <- GenomeAxisTrack()

plotTracks(list(idxTrack,axisTrack,Bam_file), from=840868,to=1002890) #,showTitle = FALSE
```
![gviz-bam](https://github.com/caokai001/Bioinformatic_code/blob/master/%E5%9B%BE%E7%89%87%E8%B5%84%E6%96%99%E5%BA%93/gviz-bam.png)

###2.画bw 文件
```
chr1 与1 的问题： https://support.bioconductor.org/p/71215/
setwd("C:\\Users\\16926\\Desktop\\CTCF")
library(rtracklayer)
library(Gviz)
library("BSgenome.Hsapiens.NCBI.GRCh38")
options(ucscChromosomeNames=FALSE)
axTrack <- GenomeAxisTrack()

## itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr17")

## levels(itrack@bandTable$chrom) <- sub("^chr", "", levels(itrack@bandTable$chrom), ignore.case=T)

## itrack@bandTable$chrom<-"17"
allChromosomeCoverage <- import.bw("PrEC_SRR513122.deeptools.bw",
                                   as="GRanges") 

accDT <- DataTrack(allChromosomeCoverage,chomosome="17",col="red",name = "chr17.peak.bw") 

plotTracks(c(axTrack,accDT), 
           from=840868,to=1002890, table = "gc5Base", 
           chromosome="17",window = 250,type = "polygon",
           fill.histogram = "blue")
```
