# WORKING WITH R

## How to make a Venn diagram from peaks files?

Using _Vennerable_, an R package available at _R-Forge_. Its installation may depends on any of the next packages.

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
(*incomplete*...)

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
gr <- filtered$paired_reads_1
df <- data.frame(chr=as.character(seqnames(gr)),starts=start(ranges(gr)),ends=end(ranges(gr)), strand=as.character(strand(gr)),id=c(gr$id),flag=c(gr$flag))
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

Using the `crond` daemon. This program should be running since the computer was reboot.

```bash
ps -ef | grep crond
```


To get access to the schedules jobs, use the `crontab` command. The `-l` option will list all actived jobs. 

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

Using the `ln` command. You can use `df` to browse your disk file system.

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


## How to work safely in the server?

Using `screen` you can have multiple virtual terminals with just one physical process. 

```bash
screen --help
```

Once you have connected to your dedicated server, `screen` allows you to create a virtual session on multiple windows.

```bash
ssh runic
screen
```

A good practice is to invoke `screen` by parsing a session name using the `-S` option. 

```bash
screen -S myscreen
```

Instead of shell commands, screen uses <kbd>Ctrl</kbd>+<kbd>a</kbd> as a signal to send commands. To get help, just use <kbd>Ctrl</kbd>+<kbd>a</kbd> and  <kbd>?</kbd> . To exit and kill the active screen use  <kbd>Ctrl</kbd>+<kbd>a</kbd>  and   <kbd>k</kbd> . If you need to close your session, detach `screen` by using  <kbd>Ctrl</kbd>+<kbd>a</kbd>  and  <kbd>d</kbd> . Later, you can list your virtual sessions and re-attach then.

```bash
ssh ms2188@tropic
screen -ls
#There is a screen on:
#	24298.myscreen		(01/08/13 11:25:00)	(Detached)
#	24254.pts-0.tropic	(01/08/13 11:26:00)	(Detached)
#1 Socket in /var/run/screen/S-ms2188.
screen -r 24298
```

> To improve your experience using `screen`, consider to customize the environment by setting up your `.screenrc` screen. For example:  [.screenrc](https://github.com/mscastillo/bash/blob/master/.screenrc).




# PROCESSING DATA


## How to choose random regions from a genome?

If the random regions should be biologically relevant (not falling in centromeres nor telomeres), use the repeat masker table which screens DNA sequences for interspersed repeats and low complexity DNA sequences.

For getting a recent version of this table, go to UCSC table menu and choose the *rmsk* table from your selected genome.![image alt text](image_0.png)

After download it, use `awk` and `shuf` to grab the columns with the genomic coordinates and randomly choose a given number of them.

```
cat repeatmasker_mm9_all_fields.tsv | awk 'BEGIN{OFS="\t";OFMT="%.f"}{print $6,($7+$8)*0.5-200,($7+$8)*0.5+200}' > repeatmasker_mm9_peaks_400bp.bed
shuf -n 1000 repeatmasker_mm9_peaks_400bp.bed > random_peaks.bed
```


## How to digest a genome?

Using `hicup_digester` to generate a restriction fragment file in bed format. This program belongs to the HiCup suite (available at http://bioinformatics.babraham.ac.uk/projects) that requires an input genome, in fasta format, the sequence pattern and the cleavage site recognised by the enzyme of interest. 

```bash
hicup_digester -1 A^AGCTT,HindIII m*.fa
```


## How to convert FASTQ format to FASTA format?

Combining `cat` and `perl` commands.

```bash
cat file_in.fastq | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){print;$i=-3}$i++;}' > file_out.fasta
```


## How to align fastq files?

Aligners: (*i*) `bwa`, (*ii*) `bowtie2` and (*iii*) `tophat`.

### bwa

First check that you have generated an index for your mapping genome. If not, you can generate it using the `index` option.

```bash
#  copy the fasta file into the folder you want to store the index
cd /usr/local/bin/bwa/index/human/hg19
cp ~/hg19.fa .
bwa index hg19.fa
```

To perform the alignment, use the `mem` option to get directly a SAM file format.  This option will use the BWA-MEM algorithm. Alternatively, you can use `bwasw` to have additional penalty and gap options. Use the `-t` option to control the number of cores to use.

```bash
bwa mem -t 4 -R "@RG\tID:mysample" /usr/local/bin/bwa/index/human/hg19/hg19.fa ~/data/mysample.fastq > mysample.sam
```


When having paired-end reads, provided the two mates of each pair in two single files.

```bash
bwa mem -t 4 -R "@RG\tID:mypesample" /usr/local/bin/bwa/index/human/hg19/hg19.fa ~/data/mypesample_1.fastq ~/data/mypesample_2.fastq  > mypesample.sam
```

### bowtie2

Once again, generate the corresponding index to your mapping genome by using `bowtie2-build`.

```bash
# copy the fasta file into the folder you want to store the index
cd /usr/local/bin/bowtie2/index/human/hg19
cp ~/hg19.fa .
bowtie2-build hg19.fa hg19.fa
```

Once your index is ready, run `bowtie2` with `-q` option for inputs in fastq format and parse the index by using the `-x` option. For multicore options use -`p`. Use `-U` to parse a coma separated list of inputs (or a single file) and `-S` to have the output in *sam* format.

```bash
bowtie2 -q -p 4 -x /usr/local/bin/bowtie2/index/human/hg19/hg19.fa -U ~/data/mysample.fastq -S ~/data/mysample.sam
```

### tophat2

Top hat uses bowtie indexes. See how to create it above.

```bash
#  copy the fasta file into the folder you want to store the index
cd /usr/local/bin/bowtie2/index/human/hg19
cp ~/hg19.fa .
bowtie2 hg19.fa hg19.fa
```


## How to split paired-end reads from a BAM/SAM file?

Using the `samtools` with the `view` option for extracting a subset of reads from a *sam* or a *bam* file. To filter out the reads we are not interested in, use the `-f` option and specify which flag (in hexadecimal format) should contain the reads to keep. For the first and second pairs, use the *0x0040* and *0x0080* flags respectively. Use the `-bh` options to get a BAM format output with its corresponding header.

```bash
samtools view -bh -f 0x0040 sample.bam > sample_paired_reads_1.bam
samtools view -bh -f 0x0080 sample.bam > sample_paired_reads_2.bam
```


## How to visualize long-range chromosomal interactions?

Using the *WashU Epigenome Browser* ([WUEB](http://epigenomegateway.wustl.edu/)). The WUEB requires paired-end bed file (tab-separated) with the next format: (*i*) chromosome name, (*ii*) start, (*iii*) end,  (*iv*) coordinate of the mate and a score of the interaction separated by commas, (*v*) identifier (a unique non-negative integer) and (*vi*) the strand (use a dot if unknown it).

```bash
echo 'chr1   111   222   chr2:333-444,55   1   .' >  mybed.bed
echo 'chr2   333   444   chr1:111-222,55   2   .' >> mybed.bed
echo 'chr3   777   888   chr2:777-888,31   1   +' >> mybed.bed
echo 'chr3   555   666   chr1:555-666,31   2   -' >> mybed.bed
```

The *bed* file should be sorted and compressed with `bgzip`. Use the `bedSort` tool from the [UCSC programs](http://hgdownload.cse.ucsc.edu/admin/exe/) instead of using any other sorting function. Then, use `tabix` to generate an index of the sorted and compressed bed file.

```bash
bedSort mybed.bed mybed.sorted.bed
cat mybed.sorted.bed
#chr1   111   222   chr2:333-444,55   1   .
#chr2   333   444   chr1:111-222,55   2   .
#chr3   555   666   chr1:555-666,31   2   -
#chr3   777   888   chr2:777-888,31   1   +
bgzip mybed.sorted.bed
tabix -p bed mybed.sorted.bed.gz
```

Later, upload both files to your server and load the data track into the WUEB by parsing the url of the compressed bed file.


## How to process ChIP-seq data from GEO?

Use `get_data` to process samples from GEO. This script is available on tropic and runic.

```bash
ssh tropic
get_data --help
```

This utility will automatically download the *sra* files and transform them to raw *fastq* format. Then, it will run `fastqc` to check the quality on each *fastq* file. It will trim the adapters and overrepresented sequences (if present). Later, it will merge the *fastq* files and will run `bowtie2`.

```bash
get_data -g GSE26014 -m GSM638307 -s SRX/SRX038/SRX038907 -x hg10
```
`get_data` will create a deep structure of folders to store the output files in your current directory.

> Before start processing any sample, check which server has less jobs running and if you have enough free disk space by using  `top` (use <kbd>q</kbd> to quit) or  `ps -aF` and `df -h`.


## How to process a batch of ChIP-Seq data?

Writing a shell script with all the command lines. The script must start with the *shebang* interpreter line. Subsequently, add the commands you want to execute. For a large batch, consider to save the standard outputs into a log file.

```bash
#! /bin/bash
date > batch.log
# first sample
get_data -g GSE26014 -m GSM638307 -s SRX/SRX038/SRX038907 -x hg19 >> batch.log
sam2bigWIG -g GSE26014 -m GSM638307 -x hg19 >> batch.log
date >> batch.log
# second sample
get_data -g GSE26014 -m GSM638308 -s SRX/SRX038/SRX038908 -x hg19 >> batch.log
date >> batch.log
sam2bigWIG -g GSE26014 -m GSM638308 -x hg19 >> batch.log
date >> batch.log
```

Before execute your script, check whether it is executable.

```bash
chmod u+x batch.sh
./batch.sh
```

Alternatively, you can run it in the background and monitor the progress by displaying the log file.

```bash
./batch.sh &
tail -f batch.log
```

> For jobs that might take a long time to finish, it is highly recommended the use of a screen.
