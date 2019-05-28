### 1.画nam 文件
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
