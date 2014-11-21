
# WORKING WITH R

## How to make a Venn diagram from peaks files?

Using _Vennerable_, an R package available at_ R-Forge_. Its installation may depends on any of the next packages.

```r
source("http://bioconductor.org/biocLite.R")
install.packages(c("graph", "RBGL"), dependencies=TRUE)
install.packages("Vennerable", repos="http://R-Forge.R-project.org")
```

After install it, load the required libraries and import your data.

```r
library(IRanges)
library(Vennerable)
library(GenomicRanges)
peaks = read.table( file="GSM0001.bed" )
colnames(peaks) = c("chr","start","end")
mygrange = GRanges( seqnames=peaks$chr,range=IRanges(start=peaks$start,end=peaks$end,names=paste(peaks$chr,peaks$start,peaks$end,sep="_")),strand="*" )
```


## How to install an older package?

You can get your current version of _R_ by using the `--version` option.


```r
R --version
#  R version 3.0.1 (2013-05-16) -- "Good Sport"
```


To manually install an old package not supported by your current version use the terminal command with the `INSTALL` option.

```r
R CMD INSTALL myoldpacakage.tar.gz
```


## How to create a _GRanges_ from a _bed_ file?

_GRanges_ (Genomic Ranges) is a class able to record genomic intervals. It allows to store basic information from a _bed_ file type (chromosome, start/end coordinates and strand) as well as additional metadata (IDs, scores, etc...). To use it, first load the package and import the _bed_ file using `read.table`. Later, give a name to the columns of your table and define the _GRanges_ with it by using the `with` function.

```r
library('GenomicRanges')
table = read.table( file='sample_1068_paired_reads_1.resampled.bed' )
colnames(table) = c('chr','start','end','strand','id','score')
paired_reads_1 = with(table,GRanges(chr, IRanges(start, end), strand, id=id, flag=score))
```

You can also stored several `GRanges` denoting paired-end reads in a single object.

```r
table = read.table( file='sample_1068_paired_reads_2.resampled.bed' )
colnames(table) = c('chr','start','end','strand','id','score')
paired_reads_2 = with(table,GRanges(chr, IRanges(start, end), strand, id=id, flag=score))
reads = list( paired_reads_1,paired_reads_2 )
names(reads) = c('paired_reads_1','paired_reads_2')
```


## How to create a bed file from a _GRanges_?

First, create your own _GRanges_ variable (find it out [here](#bookmark=id.3s34tsw03pim)) or load it from any _.Rda_ file using `load`.

```r
load('genomic_ranges.Rda')
ls()
#  [1]   "paired_reads_1"   "paired_reads_2"
class(paired_reads_1)
#  [1] "GRanges"
#  attr(,"package")
#  [1] "GenomicRanges"
paired_reads_1[1]
#  GRanges with 1 range and 2 metadata columns:
#        seqnames               ranges strand |                                         id      flag
#           <Rle>            <IRanges>  <Rle> |                                <character> <integer>
#    [1]       14 [73565650, 73565700]      - |    KL136:231:D28EKACXX:6:1101:10000:100839       115
```


Create a _data-frame_ with your selected features-metadata and write it into a file using `write.table`.

```r
gr = filtered$paired_reads_1
df = data.frame(chr=as.character(seqnames(gr)),starts=start(ranges(gr)),ends=end(ranges(gr)), strand=as.character(strand(gr)),id=c(gr$id),flag=c(gr$flag))
write.table(df,file="paired_reads_1.bed",quote=F,sep="\t",row.names=F,col.names=F)
```


# WORKING ON THE SERVERS

## How to mount a server as a local drive?

Using the `ssh` file system client. In the next example, the folder called _codex_ is mounted as a local drive under the name _superhanz_.

```bash
sshfs -o idmap=user ms2188@superhanz.cscr.cam.ac.uk:codex/ superhanz
```

## To limit resources to any feature, use the corresponding flag followed by the limit you want to set. For instance, if you want to limit the memory by the half of the maximum available, use `-m` followed by your limit in kilobytes. You can check the total and free memory by using the `free` command. Use `unlimited` to restore it to default.

```bash
free
#               total       used       free     shared    buffers     cached
#  Mem:      60000000   10000000   50000000          0      19900    9520784
#  -/+ buffers/cache:   42010312   23954420
#  Swap:     78977020      53180   78923840
ulimit -m 30000000 ; ./script.sh ; ulimit -m unlimited
```

