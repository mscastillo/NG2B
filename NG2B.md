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

```bash
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


## How to limit the amount of memory to use?

Using the `ulimi` command you can control the resources of the server. This command is available on package `util-serve`. Use `-a` option to list all possible features to be set. To get the current value of any feature, just use the option without parsing any value.

```bash
sudo aptitude search util-vserver
ulimit -a
#  core file size          (blocks, -c) 0
#  data seg size           (kbytes, -d) unlimited
#  …
ulimit -d
#  unlimited
```

To limit resources to any feature, use the corresponding flag followed by the limit you want to set. For instance, if you want to limit the memory by the half of the maximum available, use `-m` followed by your limit in kilobytes. You can check the total and free memory by using the `free` command. Use `unlimited` to restore it to default.

```bash
free
#               total       used       free     shared    buffers     cached
#  Mem:      60000000   10000000   50000000          0      19900    9520784
#  -/+ buffers/cache:   42010312   23954420
#  Swap:     78977020      53180   78923840
ulimit -m 30000000 ; ./script.sh ; ulimit -m unlimited
```


## How to set scheduled jobs?

Using the _crond_ daemon. This program should be running since the computer was reboot.

```bash
ps -ef | grep crond
```


To get access to the schedules jobs, use the _crontab_ command. The _-l_ option will list all actived jobs. 

```bash
crontab -l
```


To edit it, use the `-e` option. Add a new line to set a new scheduled job. The timestamp format is:

## How to backup a MySQl database?

Using the `mysqldump` command.

```bash
mysqldump -u $USER -p$PASSWORD $STAGINGDB > staging.sql
mysqldump -u $USER -p$PASSWORD $BIOINFORMATICSDB > bioinformatics.sql
```


## How to run parallel jobs?

Using the `parallel` command. Use brackets in the statement as the input source, with the arguments placed after the triple colon.

```bash
parallel -help
parallel 'gunzip {}' ::: ls *-*/*.fq.gz
parallel 'fastqc --quiet -f fastq {}' ::: ls *-*/*.fastq
parallel 'gzip -f {}' ::: ls */*-*.fastq
```


## How to create a symbolic link?

Using the **ln** command. You can use **df** to browse your disk file system.

```bash
df -h
#Filesystem            Size  Used Avail Use% Mounted on
#/dev/sdc1              37G  5.5G   30G  16% /
#/dev/sdb1             2.7T  1.9T  672G  75% /projects
#/dev/sda1             2.7T  2.4T  234G  92% /home
#//pacific/gottgens/    22T   21T  963G  96% /home/rlh60/Pacific
#//pacific/huntly/      19T  9.1T  9.2T  50% /home/rlh60/Brian_Pacific
#apollo:/export/data    19T  5.7T   12T  33% /data
#...
cd
ln -s /data
ls -l
```


You can also create symbolic links of your favourites programs into **/bin** to avoid using the full path.

```bash
# do not use the full path to your programs
/home/Programs/bowtie-0.12.9/bowtie -m 1 -v 2 -S --phred33-quals hg19_ucsc A006.fastq > A006.sam
# use a symbolic link instead
ln -s /home/Programs/bowtie-0.12.9/bowtie /bin
bowtie -m 1 -v 2 -S --phred33-quals hg19_ucsc A006.fastq > A006.sam
```


## How to manage RSA keys for authentication?

Using `ssh-keygen` you can generate private-public key pairs. Don’t use a _passphrase_ to avoid user monitoring every time you want to login.

```bash
ssh ms2188@runic
ssh-keygen -t rsa -C "RSA key from ms2188@runic"
#Created directory '/home/ms2188/.ssh'.
#Enter passphrase (empty for no passphrase):
#Enter same passphrase again:
#Your identification has been saved in /home/ms2188/.ssh/id_rsa.
#Your public key has been saved in /home/ms2188/.ssh/id_rsa.pub.
#...
cat .ssh/id_rsa.pub
```

Now just add your _public key_ into the _authorized keys_ file in the server of interest.

```bash
ssh ms2188@superhanz
nano .ssh/authorized_keys
```

> To speed up the connection between servers, consider to set an alias such as `alias superhanz="ssh ms2188@superhanz"`.
