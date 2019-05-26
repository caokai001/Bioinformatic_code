library(GenomicRanges)
library(rtracklayer)
library(devtools)
library(derfinderData)
bw_file="SRR1282236.trim.deeptools.bw"
gr=GRanges(c("1","21"),IRanges(c(9481564,9963234),width = 1000))
A=import(BigWigFile(bw_file),selection=gr[1])

###https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/viz_course/Presentations/singlepage/Viz_part_2.html#/vizdata
###加载
library(rtracklayer) 
allChromosomeCoverage <- import.bw("SRR1282236.trim.deeptools.bw",
                                   as="GRanges") 
options(ucscChromosomeNames=FALSE)
accDT <- DataTrack(allChromosomeCoverage,chomosome="17") 
plotTracks(c(accDT), 
           from=1000000,to=16000000, 
           chromosome="17",type="hist") 
