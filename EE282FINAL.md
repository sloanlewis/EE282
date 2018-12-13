# EE282 FINAL PROJECT
## Background
The premise of this project is to discover the effects of chronic alcohol consumption on the mucosal immune system in the gut. We utilize a nonhuman primate model of voluntary ethanol self-administration where rhesus macaques have open access to a 4% w/v ethanol solution for 12 months. After these 12 months, the monkeys are necropsied and biopsies from the gut are taken. We isolated mononuclear cells from the lamina propria of the duodenum, jejunum, ileum, and colon of control (n=4) and ethanol-drinking (n=8) rhesus macaques. Once these cells were isolated, they were plated and stimulated with PMA/ionomycin for 16 hours. The supernatants were taken for luminex analysis and the cells for RNA-seq. We extracted RNA, built libraries and then multiplexed libraries were subjected to single-end 86 or 100 base pair sequencing using the Illumina NextSeq 500 (jejunal and ileal LPL) or HiSeq 4000 (duodenal and colonic LPL) platforms, respectively. This is the data I will be using for this project.

Note: In the Messaoudi Lab, we use the UC Riverside HPCC so I will try to make everything run on the UCI cluster, but errors may be due to this.

## Set-up
```
qrsh -q abio,free128,free88i,free72i,free32i,free64 -pe openmp 32
#Create working directories for the final project
mkdir /pub/jje/ee282/sloan/Final_Project
cd /pub/jje/ee282/sloan/Final_Project
mkdir data results sequences

```

## Get the data
Because it would be too much to transfer all 68 sequence data files, I will make symbolic links to 2 of the raw sequences in the sequences directory.

```
cd sequences
ln -s ~/Final_Sequences/*.gz ./
#to see the sequences
ll

25787_stim_ile_LPL.read1.fastq.gz -> /data/users/sloanal/Final_Sequences/25787_stim_ile_LPL.read1.fastq.gz
25787_un_ile_LPL.read1.fastq.gz -> /data/users/sloanal/Final_Sequences/25787_un_ile_LPL.read1.fastq.gz

```
## Sequence quality check
Run fastqc on the sequences to check their quality and see how much should be trimmed.

```
#run fastqc on both samples, if you have a lot make a script to run on a whole directory
module load jje/jjeutils
module load jje/kent
module load fastqc
fastqc *.fastq.gz

```
Check the fastqc report and make sure the sequences are of adequate quality. This website explains each parameter: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Use this to determine how much needs to be trimmed off of each sequence.


## Trimming
Sequences need to be trimmed to get rid of adapters added for library preparation and multiplexing.

```
module load trimmomatic/0.35
#Set parameters and run on each file
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 file_name.fastq.gz file_name_trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:4 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#re-run fastqc after trimming to make sure quality is still good and the adapters were trimmed properly
fastqc *trimmed.fq.gz

```

## Alignment to the macaca mulatta genome

__Downloading:__

```
#Download the most recent version of the macaca mulatta genome from ensembl along with the GTF and CHECKSUMS file
#genome
wget ftp://ftp.ensembl.org/pub/release-94/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa.gz
#genome CHECKSUMS
wget ftp://ftp.ensembl.org/pub/release-94/fasta/macaca_mulatta/dna/CHECKSUMS
mv CHECKSUMS CHECKSUMS_genome
#GTF annotation file
wget ftp://ftp.ensembl.org/pub/release-94/gtf/macaca_mulatta/Macaca_mulatta.Mmul_8.0.1.94.gtf.gz
#GTF CHECKSUMS
wget ftp://ftp.ensembl.org/pub/release-94/gtf/macaca_mulatta/CHECKSUMS

```

__Quality of download check:__

```
sum Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa.gz
#compare to CHECKSUMS file

sum Macaca_mulatta.Mmul_8.0.1.94.gtf.gz
#compare to CHECKSUMS.1 file

```

I then ran the alignment using Dr. Girke's systemPipeR pipeline: https://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeRNAseq.pdf
I will not go through the whole thing for this project as it is already published. After aligning the sequenced reads to the genome using Tophat, the mapped reads are counted using the annotation file. 

## Data Plots
I put my final counts and rpkm files as well as the targets file into the results folder to use for the following plots.

__Principle Component Analysis: 2 ways__

*To look at only the ileum samples (controls, moderate and heavy drinkers, stim, and nostim):*

```
cd ../results
nano targets_gut.txt
#put a # in front of the lines you don't want to use, in this case everything except the ileum samples
module load R
R
library(DESeq2)
#read in files
countDF <- read.delim("./countDF_MASTER.xls", row.names=1)
targets <- read.delim("targets_gut.txt", comment.char = "#")
#subset counts file to only use the samples from the targets file
names <- as.character(targets$SampleName)
countDF2 <- subset(countDF, select=names)
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF2, colData = colData, design = ~ condition)
rld <- rlog(dds)
pdf("./PCA_group_rlog_ileum.pdf")
plotPCA(rld)
dev.off()
q()

```

MC_INS= male control ileum nostim

MM_INS= male moderate drinker ileum nostim

MH_INS= male heavy drinker ileum nostim

![PCA_ileum](https://github.com/sloanlewis/EE282/blob/master/PCA_group_rlog_ileum.pdf)

*To look at all gut sections only nostim samples to see if there are any innate differences in transcriptional profiles of lamina propria lymphocytes from different sections, controls and drinkers combined:*

```
nano targets_gut.txt
#put a # in front of the lines you don't want to use, in this case all stim samples
library(DESeq2)
countDF <- read.delim("./countDF_MASTER.xls", row.names=1)
targets <- read.delim("targets_gut.txt", comment.char = "#")
#subset counts file to only use the samples from the targets file
names <- as.character(targets$SampleName)
countDF2 <- subset(countDF, select=names)
colData <- data.frame(row.names=targets$SampleName, condition=targets$SectionOnly)
dds <- DESeqDataSetFromMatrix(countData = countDF2, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("./PCA_group_rlog_bysection.pdf")
plotPCA(vsd)
dev.off()

```
![PCA_group](https://github.com/sloanlewis/EE282/blob/master/PCA_group_rlog_bysection.pdf)

*For the purposes of this project, I am setting a high cutoff on the counts file to make a heatmap with fewer genes and without running DEG analysis first. I am looking for expression differences in these highly expressed genes in the drinkers versus controls in the ileum(nostims)*

```
module load R
R
rpkm <- read.delim("./rpkmDFeByg_MASTER.xls", row.names=1)
# grab only ileum nostims
r <- rpkm[,c(20, 21, 22, 23, 24, 25, 26, 27, 28)]
#re-order for the heatmap
s <- r[, c(6,7,8,1,3,5,2,4,9)]
u <- cbind(rownames(s), data.frame(s, row.names=NULL))
#set rowmeans cutoff for counts file
counts <- read.delim("./countDF_MASTER.xls", row.names=1)
c <- subset(counts, rowMeans(counts)>=100)
d <- cbind(rownames(c), data.frame(c, row.names=NULL))
f <- subset(d, select=c(1))
# merge the genes from the counts file with the cutoff with your rpkm file subsetted to have only the nostim ileum samples
merge <- merge(f, u, by.x = "rownames(c)", by.y = "rownames(s)")
final <- data.frame(merge[,-1], row.names=merge[,1])

#Run this for a clustered heatmap
args <- commandArgs(TRUE)
library("RColorBrewer")
library(gplots)
library(lattice)

#Define the colors
hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(100)

# file1 is the list of merged genes
file1  <- args[1]
# Output name - heatmap pdf prefix
outprefix <- args[2]

#Open file handles and save the first columns
f1 <- final

############################################################
y <- t(scale(t(as.matrix(f1))))

#Clustered
pdf(paste(outprefix, "_clustered_heatmap.pdf", sep=""))
heatmap.2(as.matrix(f1), col=rev(hmcol), scale="row", key=T, keysize=1.5, density.info="density", cexRow=0.5, cexCol=0.9, Colv=FALSE, distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="average"), trace='none')
dev.off()

```
![heatmap](https://github.com/sloanlewis/EE282/blob/master/NA_clustered_heatmap.pdf)
