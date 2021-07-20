# Daphnia Resting Eggs Lake Constance

## Pre-processing sequencing data

### Eva

Put in here what you have done from raw data to mapped reads.
Include tools, versions, parameters, etc...










### Tania

Note: From sam files to read depth we used the Snakemake workflow from the [RestEggs_snakemake_pbs Tutorial](https://github.com/tholtzem/RestEggs_snakemake_pbs). The analysis is based on the [physalia-lcwgs: the Physalia course on Population genomic inference from low-coverage whole-genome sequencing data, Oct 19-22 2020](https://github.com/nt246/physalia-lcwgs).

The snakemake workflow was performed on mach2 using Anaconda3/2020.03/miniconda-base-2020.03 and snakemake-minimal=6.


Snakemake rules used in this workflow:

1) sam2bam

tool: samtools view
version: samtools=1.12

```
rule sam2bam:
  input:
    sam = '/home/uibk/c7701125/scratch/EXCHANGE/RESTEGGS_LC/MAPPING/bbmap/{sample}.sam.gz'
  output:
    bam = 'bbmap/rapid/{sample}.bam'
  message: """--- Converting sam files to sorted bam ---"""
  threads: 12 
  shell:
    """
    samtools view -F 4 -Shu {input.sam} | samtools sort - -o {output.bam}
    """
```

2) Then I included data from our re-sequencing data set.

I selected 9 clonal lines for *D. longispina* , *D. galeata*, *D. cucullata* (3 clones for each species).

```
rsync -avP --files-from=ref_clones.txt /home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/bbmap/rapid/ bbmap/rapid/

```

This data was already pre-processed including mapping and merging of mapped reads, see ref_clones.txt for a list of these files. For more details on the pre-processing steps of the re-sequencing data set, see the Snakemake workflow from the [daphnia_snakemake_pbs Tutorial](https://github.com/tholtzem/daphnia_snakemake_pbs).

3) Filter reads with mapping quality lower than 20

tool: samtools view
version: samtools=1.12

```
rule bam_q20:
  input:
    bam = 'bbmap/rapid/{sample}.bam'
  output:
    bamq20 = 'bbmap/rapid/minq20/{sample}.minq20.bam'
  message: """--- Filter bam reads with a mapping quality lower than 20 and convert to sorted bam ---"""
  threads: 12 
  shell:
    """
    samtools view -Shu -q 20 {input.bam} | samtools sort - -o {output.bamq20}
    """
```

4) Remove duplicates

tool: picard MarkDuplicates
version: picard=2.25.0-0

```
rule remove_duplicates:
  input:
    bam = 'bbmap/rapid/minq20/{sample}.bam'
  output:
    deDup = 'deDup/{sample}.dedup.bam',
    metrics = 'deDup/{sample}.dedup.metrics.txt'
  log: 'log/{sample}.dedup.bam.log'
  threads: 12
  message: """--- Removing duplicates of merged bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.bam} OUTPUT={output.deDup} METRICS_FILE={output.metrics} 2> {log}
    """
```

5) Mark duplicates

```
rule mark_duplicates:
  input:
    mrgd = 'bbmap/rapid/minq20/{sample}.bam'
  output:
    mrkDup = 'markDup/{sample}.mrkdup.bam',
    metrics = 'markDup/{sample}.mrkdup.metrics.txt'
  log: 'log/{sample}.mrkdup.bam.log'
  threads: 12
  message:
    """--- Marking duplicates of merged bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.mrkDup} METRICS_FILE={output.metrics} 2> {log}
    """
```

6) Clip overlaps

tool: bam clipOverlap
version: bamutil=1.0.15-0

```
rule clip_overlap:
  input:
    deDup = 'deDup/{sample}.dedup.bam'
  output:
    clip = 'deDup/{sample}.overlapclipped.bam' 
  log: 'log/{sample}.overlapclipped.bam.log'
  threads: 12
  message:
    """ Clip overlapping paired end reads """
  shell:
    """
    bam clipOverlap --in {input.deDup} --out {output.clip} --stats 2> {log}
    """
```

7) Index clipped bam files (required by gatk)

tool: samtools index
version: samtools=1.12

```
rule Index_clippedBAM:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam'
  output:
    idx = 'deDup/{sample}.overlapclipped.bam.bai'
  log: 'log/{sample}.overlapclipped.bam.log'
  threads: 2
  message: """--- Indexing clipped BAM files with samtools ---"""
  shell:
    """
    samtools index {input.clip} {output.idx} 2> {log}
    """
```

8) Index reference genome (required by gatk)

tool: samtools faidx
version: samtools=1.12
reference genome: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta

```
rule refIndex:
  input:
    ref = config['ref']
  output:
    'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai'
  log: 'log/refIndex.log'
  shell:
    """
    samtools faidx {input.ref} 2> {log}
    """
```

9) Create dictionary of reference genome (required by gatk)

tool: picard=2.25.0-0
picard CreateSequenceDictionary
reference genome: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta

```
rule ref_Dict:
  input:
    ref = config['ref']
  output:
    'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dict'
  log: 'log/refDict.log'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar CreateSequenceDictionary R={input.ref} O={output} 2> {log}
    """
```

10) Create list of clipped bam files

```
rule ls_ClipBam:
  input:
    #bam = 'bbmap/all/{sample}.bam'
  output:
    #touch('bbmap/all/{sample}.ls.done'),
    'list/overlapclippedBAM.list'
  log: 'log/lsBam.log'
  message: """--- Creating a sample list of clipped bam files for indel-realignement ---"""
  shell:
    """
    ls deDup/*minq20.overlapclipped.bam > {output} 2> {log}
    """
```

11) Create list of potential indels across all samples

tool: gatk RealignerTargetCreator
version: gatk=3.8

requirements:
* reference genome: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta
* indexed reference: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai
* dictionary of reference: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dict

```
rule list_indels:
  input:
    #clip = 'deDup/{sample}.overlapclipped.bam',
    clip = 'list/overlapclippedBAM.list',
    ref = config["ref"]
  output:
    #intervals = 'list/{sample}.intervals'
    intervals = 'list/all_samples.intervals'
  log: 'log/listIndels.log'
  threads: 12
  message:
    """ Create list of potential in-dels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx120g -jar $GATK -T RealignerTargetCreator -R {input.ref} -I {input.clip} -o {output.intervals} -drf BadMate 2> {log}
    """
```
 
12) indel-realignment

tool: gatk IndelRealigner
version: gatk=3.8

requirements:
* reference genome: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta
* indexed reference: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai
* dictionary of reference: ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dic

```
rule realign_indel:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam',
    ref = config["ref"],
    intervals = 'list/all_samples.intervals'
  output:
    realigned = 'realigned/{sample}.realigned.bam'
  log: 'log/{sample}.realigned.bam.log'
  threads: 12
  message:
    """ Realign in-dels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx120g -jar $GATK -T IndelRealigner -R {input.ref} -I {input.clip} -targetIntervals {input.intervals} -o {output.realigned} --consensusDeterminationModel USE_READS 2> {log}
    """
```

13) Estimate read depth
tool: samtools depth
version: samtools=1.12

```
rule bam_depth:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    depth = 'depth/{sample}.realigned.bam.depth.gz'
    
  threads: 12
  message:
    """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input.realigned} | cut -f3 | gzip > {output.depth}
    """
```

### Read depth stats

14) List samtools depth files

```
rule ls_depth:
  output:
    'list/depth.list'
  log: 'log/depth.list.log'
  message: """--- Creating a sample list of samtools-depth files for Rscript ---"""
  shell:
    """
    ls depth/*depth.gz | cut -f2 -d '/' > {output} 2> {log}
    """
```

15) Read depth stats ALL samples (to get an overview of our data)

* summary stats table for ALL samples
* plot overall depth distribution and depth summary for chromosom "dgal1" only  

```
rule read_depth:
  #input:
    #test.list
  output:
    #pdf = "some.pdf"
    touch('genome_stats.done')
  log: 'log/genome_stats.done.log'
  threads: 12
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/read_depth_R.R 2> {log}
    """
```

Rscript:

```

library(tidyverse)

getwd()
setwd("/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs/")

basedir <- "depth/" # Make sure to edit this to match your $BASEDIR
bam_list <- read_lines(paste0("list/depth.list"))

#bam_list

getwd()
for (i in 1:69){
    bamfile = bam_list[i]
    # Compute depth stats
    depth <- read_tsv(paste0(basedir, bamfile), col_names = F)$X1
    mean_depth <- mean(depth)
    sd_depth <- sd(depth)
    mean_depth_nonzero <- mean(depth[depth > 0])
    mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
    median <- median(depth)
    presence <- as.logical(depth)
    proportion_of_reference_covered <- mean(presence)
      
  # Bind stats into dataframe and store sample-specific per base depth and presence data
  if (i==1){
    output <- data.frame(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered)
    total_depth <- depth
    total_presence <- presence
  } else {
    output <- rbind(output, cbind(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered))
    total_depth <- total_depth + depth
    total_presence <- total_presence + presence
  }
}

output %>%
  mutate(across(where(is.numeric), round, 3))
                
write.table(output,"depth_statistics.txt", sep ="\t", quote = F)

#set genome coordinates for plots, the depth files are a list of depth values for each site. Genomic coordinates for scaffolds can be e.g. calculated from the *.fai file (first column = scavvold name, second column = scaffold length). For example: gal1= 1:2950711; gal2 = 2950712-5869635...) 
#dgal1
chr <- "dgal1"
coord_start <- 1
coord_end <- 2950711
total_depth <- total_depth[coord_start:coord_end]
total_presence <- total_presence[coord_start:coord_end]

#Plot the depth distribution

pdf(file=paste0(basedir,chr,"_depth_distribution.pdf"))
#tibble(total_depth = total_depth, position = 1:length(total_depth))  %>%
tibble(total_depth = total_depth, position = coord_start:coord_end)  %>%
  ggplot(aes(x = position, y = total_depth)) +
  geom_point(size = 0.1)
dev.off()

# Total depth per site across all individuals 
total_depth_summary <- count(tibble(total_depth = total_depth), total_depth)
total_presence_summary <- count(tibble(total_presence = total_presence), total_presence)

pdf(file=paste0(basedir,chr,"_depth_summary.pdf"))
total_depth_summary %>%
  ggplot(aes(x = log(total_depth), y = n)) +
  geom_point()
dev.off()
pdf(file=paste0(basedir,chr,"_presence_summary.pdf"))
total_presence_summary %>%
  ggplot(aes(x = total_presence, y = n)) +
  geom_col()
dev.off()

```

16) Plot read depth using a summary table

* summary table: 'depth/depth_statistics.txt'
* A bar plot of mean read depth per sample was created to get a general idea of our data
* dataframe with samples with a mean_depth > 10 were subsampled to create a new dataframe and bar plot
* Samples having a mean read depth higher than 10 were listed in a text file

Rscript:

```
library(tidyverse)

getwd()
setwd("/mnt/data/snakemake_mach2/Daphnia_RestEggs_snakemake_pbs/")

stats <- read.table(file = 'depth/depth_statistics.txt', sep = ',', header = TRUE)

# subset dataframe stats
df <- subset(stats, mean_depth > 1)
df

dfu1 <- subset(stats, mean_depth < 1)
dfu1$bamfile

df10 <- subset(df, mean_depth > 10)
df10$bamfile

df20 <- subset(df, mean_depth > 20)
df20$bamfile

df10_u20 <- subset(df10, mean_depth < 20)
df10_u20$bamfile

# Bar lot of mean read depth per sample
#pdf("depth/stats/depth_hist_perSample.pdf")
pdf("depth/plots/depth_df10_hist_badSampleOut.pdf")
barplot(df10$mean_depth, names=df10$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth", line = 3)
#title("Mean read depth per sample")
title("Mean read depth per sample (df10)")
dev.off()

# Check distributions of all samples using boxplots
#Rs boxplot function uses the standard rule to indicate an observation as a potential outlier
#if it falls more than 1.5 times the IQR (Inter-Quartile Range, calculated as Q3-Q1) below Q1 or above Q3.
#The potential outliers are plotted with circles and the Whiskers (lines that extend from Q1 and Q3 typically to the minimum and maximum) are shortened to only go as far as observations that are within 1.5*IQR of the upper and lower quartiles. 
#The box part of the boxplot is a box that goes from Q1 to Q3 and the median is displayed as a line somewhere inside the box
pdf("depth/plots/depth_df10_boxplot.pdf")
boxplot(df10$mean_depth)
title("Distribution of mean read depth (df10)")
dev.off()

df10_sum <- summary(df10$mean_depth)
df10_sum

meanDepth_allInd <- mean(df10$mean_depth)
medianDepth_allInd <- median(df10$mean_depth)

IQRT<-IQR(df10$mean_depth)
IQRT
#IQRT<-(26.92935-19.27842)
#IQRT

# number of individuals
N = nrow(df10)
N

#1.5 times the interquartile range from the mean (median) were excluded as they are expected to be enriched for paralogs
## mean
MaxDepth1 = (meanDepth_allInd + IQRT*1.5)*N
MinDepth1  = (meanDepth_allInd - IQRT*1.5)*N
## mediam
MaxDepth2 = (medianDepth_allInd + IQRT*1.5)*N
MinDepth2  = (medianDepth_allInd - IQRT*1.5)*N
# d + 3*sqrt(d)
HengLi1_max <- (meanDepth_allInd + 3*sqrt(meanDepth_allInd))*N
HengLi2_max <- (medianDepth_allInd + 3*sqrt(meanDepth_allInd))*N

HengLi1_min <- (meanDepth_allInd - 3*sqrt(meanDepth_allInd))*N
HengLi2_min <- (medianDepth_allInd - 3*sqrt(meanDepth_allInd))*N

# alternative approach

## sum mean_depth for all individuals
total <- mean(df10$mean_depth)
# sum the square root (SQRT) of the mean depth
SQRT <- sqrt(df10$mean_depth)
# sum IQRT for all individuals
iqr <- IQR(df10$mean_depth)

## Filter1
maxFilter1 <- (total + 3*SQRT)
minFilter1 <- (total - 3*SQRT)
## Filter2
maxFilter2 <- (total + iqr*1.5)
minFilter2 <- (total - iqr*1.5)



############################################################
#non-zero

dfNonZero <- summary(df10$mean_depth_nonzero)
dfNonZero

# Bar lot of mean read depth per sample
## samples having zero depth excluded
pdf("depth/plots/depth_NonZero_hist_perSample.pdf")
barplot(df10$mean_depth_nonzero, names=df10$id,las=2,cex.names=0.40, ann = F)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Mean read depth (df10)", line = 3)
title("Mean read depth (non-zero) per sample (df10)")
dev.off()

##########################################################
# write samples having a mean read depth higher than 10 to a text file, without using quotes
#realigned <- paste0(df10$bamfile,".minq20.realigned.bam")
#Extract Characters Before Pattern using sub()
## first pattern, then replacement of the pattern, then string
realignedBAM_df10 <-sub(".depth.gz", "", df10$bamfile)
write.table(realignedBAM_df10, file = "realignedBAM_depth10.list", sep = "\n",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

```

Subselect samples with a mean read depth equal or over 10, and write respective samples to a list, named "realignedBAM_depth10.list".

Use this list as input in ANGSD.



### Popualtion structure using unpruned GLs, UPDATED 09.07.2021

#### ANGSD: Genotype likelihoods and filtering

Calculate genotype likelihoods (beagle) using various filters, for bam files having a read depth gretaer than 10.

The genotype likelihood output for these variant sites is in beagle likelihood file format (beagle.gz) by specifying the -doGlf 2 option. This will be used as input for estimating the covariance matrix using PCAngsd.

input: realignedBAM_depth10.list
output: angsd_LC_GL2_cutoff_maf0018_nInd55.mafs.gz, angsd_LC_GL2_cutoff_maf0018_nInd55.beagle.gz


```
rule angsd_GL2:
  input:
    ref = config["ref"]
  output:
    #touch('angsd/GL2_cutoffs_{IND}.done')
    touch('angsd/angsd_GL2_cutoffs_maf0018_nInd55.done')
  log: 'log/angsd_GL2_cutoffs_maf0018_nInd55.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods for bam files having a read depth gretaer than 10 and the cutoff -nInd 55 and -minMaf 0.018 (1.8 %) while recording how many sites are predicted to be variable for each scenario, """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -ref {input.ref} -out angsd/angsd_LC_GL2_cutoff_maf0018_nInd55 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 20 -minQ 20 -minInd 55 -setMinDepth 400 -setMaxDepth 2000 -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.018 -doCounts 1 2> {log}
    """
```

parameters:

-b /-bam [filelist]

-GL 2: GLs, estimated with GATK model
-doGlf 2: output in beagle format, *glf.gz
-doMajorMinor 1: Infer major and minor from GL
-doMaf 1: Estimate allele frequencies (fixed major and minor)

Note: We may be interested in looking at allele frequencies for variable sites only. To call SNPs in ANGSD we use:
-SNP_pval: Remove sites with a pvalue larger
-minMaf 0.018: Filter SNPs with a minor allel frequency, calculated using the formula 2/2*N_IND (here: 2/2*55=0.018)

-minMapQ 20: mapping quality filter
-minQ 20: base score quality filter

-setMinDepth 400: Default -1. If the total depth is below this value, the site is discarded.
-setMaxDepth 2000: Default -1. If the total depth is above this value, the site is discarded. 
-doCounts 1: Count the number A,C,G,T. All sites, All samples
Sometimes we want or need the frequency of the different bases.

-uniqueOnly 1: Remove reads that have multiple best hits (default 0).
-remove_bads 1: Remove reads tagged as bad, i.e. not primary, failure and duplicate reads (default 1).
-only_proper_pairs 1: Considering only proper pairs. 1: include only proper (default), 0: use all reads.
-trim 0: Trim [int] bases at both ends of the reads. Useful for ancient DNA (default 0).
-C 50: Reduces the effect of reads with excessive mismatches

-nInd: use only sites with data from at least N individuals (when -nInd is NOT set, -nInd 0!)


This filter is very useful, but takes a lot of time!
##-baq 1: Computes base alignment quality used to rule out false SNPs close to INDELS


output:

* .mafs.gz (allele frequencies computed for each site).
The columns are: chromosome, position, major allele, minor allele, reference allele, allele frequency, p-value for SNP calling (if -SNP-pval was called), number of individuals with data. The last column gives the number of samples with data (you can see that this never below 5 given our filtering)
Note:  allele frequency (knownEM) refers to the allele labelled as minor (a low allele frequency could reflect the fact that that site is monomorphic)

general guidance:

(https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/04_snp.md)

As a general guidance, -GL 1, -doMaf 1/2 and -doMajorMinor 1 should be the preferred choice when data uncertainty is high. If interested in analysing very low frequency SNPs, then -doMaf 2 should be selected. When accurate information on reference sequence or outgroup are available, one can use -doMajorMinor to 4 or 5. Also, detecting variable sites based on their probability of being SNPs is generally a better choice than defining a threshold on the allele frequency. However, various cutoffs and a dedicated filtering should be perform to assess robustness of your called SNPs.


Count the number of sites for default parameters

```
zcat angsd_LC_GL2_cutoff_maf0018_nInd55.mafs.gz | tail -n+2 | wc -l
```
6,430,625 sites

#### PCAngsd: Calculate a covariance matrix on unpruned data

PCAngsd can estimate a covariance matrix and individual admixture proportions simultaneously. PCAngsd can also automatically infer the most likely number of ancestry cluster or the number of clusters can be set using the -admix_K option.

```
rule PCAngsd2:
  input: 
    touched = 'angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'
  output:
    touch('pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55.done')
  log: 'log/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.beagle.gz -o pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat 2> {log}
    """
```

Plot covariance matrix
```
rule plot_covMat:
  input:
    #touched = 'pcangsd/PCAngsd_GL2_{IND}.done',
    #covMat ='pcangsd/PCAngsd_GL2_{IND}_covmat.cov'
    touched = 'pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55.done' 
  output:
    pdf = 'pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.pdf'
  log: 'log/PCAngsd_GL2_cutoffs_maf0018_plotCovMat_nInd55.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    Rscript scripts/plot_covMat.R pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.cov {output.pdf} 2> {log}
    """
```
scripts/plot_covMat.R:

```

library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
library(viridis)

getwd()
setwd("/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs/")

args <- commandArgs(trailingOnly=TRUE)

#Load the covariance matrix
cov <- as.matrix(read.table(args[1], header = F))

#We will also add a column with population assingments
sample_data <- read.csv("list/metadata_RESTEGGS_LC_OLDvsNEW_badOUT.csv", header = T, sep='\t')
sample_data$sample_id


# IDs
pop <- as.vector(sample_data$sample_id)
pop
mme.pca <- eigen(cov) #perform the pca using the eigen function
eigenvectors = mme.pca$vectors #extract eigenvectors
pca.vectors = as.data.frame(cbind(pop, eigenvectors)) #combine with our population assignments
df = type.convert(pca.vectors)
head(df)


# groups(species)

groups <- as.vector(sample_data$species)
groups
mme.pca <- eigen(cov) #perform the pca using the eigen function
eigenvectors = mme.pca$vectors #extract eigenvectors

pca.vectors2 = as.data.frame(cbind(groups, eigenvectors)) #combine with our population assignments
df2 = type.convert(pca.vectors2)
head(df2)

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
#percentage <- round((mme.pca$values/pca.eigenval.sum)*100, 2)
#percentage <- paste( colnames(df), "(", paste( as.character(percentage), "%", ")", sep="") )
varPC1 <- round((mme.pca$values[1]/pca.eigenval.sum)*100, 2) #Variance explained by PC1
varPC2 <- round((mme.pca$values[2]/pca.eigenval.sum)*100, 2) #Variance explained by PC2
varPC3 <- round((mme.pca$values[3]/pca.eigenval.sum)*100, 2) #Variance explained by PC3
varPC4 <- round((mme.pca$values[4]/pca.eigenval.sum)*100, 2) #Variance explained by PC4

PC1 <- paste( "PC1", "(", paste( as.character(varPC1), "%", ")", sep=" ") )
PC2 <- paste( "PC2", "(", paste( as.character(varPC2), "%", ")", sep=" ") )
PC3 <- paste( "PC3", "(", paste( as.character(varPC3), "%", ")", sep=" ") )
PC4 <- paste( "PC4", "(", paste( as.character(varPC4), "%", ")", sep=" ") )

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

# Colour palette with black:
cbp2 <- c("#009E73", "#000000", "#D55E00", "#56B4E9", "#E69F00",
          "#F0E442", "#0072B2", "#CC79A7")


# Open a pdf file
pdf(args[2])

## plot PCA with IDs
pca = ggplot(data = df, aes(x=V2, y=V3, label=pop)) + geom_point(pch=16, size=4.5) + ggtitle("PC1 vs PC2 (allSNPs)") + geom_text()
pca = pca + xlab(PC1) + ylab(PC2) + theme + scale_color_manual(values = cbp2)
pca

## PC1vsPC2
pca = ggplot(data = df2, aes(x=V2, y=V3, color=groups)) + geom_point(pch=16, size=4.5) + ggtitle("PC1 vs PC2 (groups, allSNPs)")
pca = pca + xlab(PC1) + ylab(PC2) + theme + scale_color_manual(values = cbp2)
pca

## PC1vsPC3
pca = ggplot(data = df2, aes(x=V2, y=V4, color=groups)) + geom_point(pch=16, size=4.5) + ggtitle("PC1 vs PC3 (groups, allSNPs)")
pca = pca + xlab(PC1) + ylab(PC3) + theme + scale_color_manual(values = cbp2)
pca

## PC2vsPC3

pca = ggplot(data = df2, aes(x=V3, y=V4, color=groups)) + geom_point(pch=16, size=4.5) + ggtitle("PC2 vs PC3 (groups, allSNPs)")
pca = pca + xlab(PC2) + ylab(PC3) + theme + scale_color_manual(values = cbp2)
pca

dev.off()

```

### Popualtion structure using LD-pruned GLs

#### ngsLD

* create geno file

```
rule prepare_genoFile:
  input:
    'angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_geno.beagle.gz'
  log: 'log/prepare_genoFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.beagle.gz | cut -f 4- | awk 'NR != 1' | awk 'NR % 50 == 0' | gzip  > {output} 2> {log}
    """
```

* create postion file

```
rule prepare_posFile:
  input:
    'angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_pos.gz'
  log: 'log/prepare_posFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.mafs.gz | cut -f 1,2 |  awk 'NR != 1' | awk 'NR % 50 == 0' | sed 's/:/_/g'| gzip > {output}
    """
```

* run ngsLD

If both --max_kb_dist and --max_snp_dist are set to 0, ngsLD will output all comparisons, even between different chromosomes/scaffolds/contigs. (https://github.com/fgvieira/ngsLD)

```
rule run_ngsLD:
  input:
   position = 'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_pos.gz',
   geno = 'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_geno.beagle.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.ld.gz'
  log: 'log/run_ngsLD.log'
  threads: 24
  message:
    """ Estimate LD using ngsLD, which employs a maximum likelihood approach to account for genotype uncertainty in low-coverage whole genome sequencing data. """
  shell:
    """
    module load ngsld/1.1.1
    ngsld --geno {input.geno} --pos {input.position} --probs --n_ind 55 --n_sites 128612 --max_kb_dist 0 --max_snp_dist 0 --rnd_sample 0.008 --n_threads {threads} | gzip --best > {output} 2> {log}
    """
```

* prune your data using prune_graph.pl

```
rule run_LDpruning:
  input:
   'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.ld.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id'
  log: 'log/LD_pruned.log'
  message:
    """ Prune your data, remove SNPs with high LD """
  shell:
    """
    module load singularity/2.x
    zcat {input} | singularity exec /apps/uibk/ngsld/1.1.1/ngsLD.sandbox /opt/ngsLD/scripts/prune_graph.pl --max_kb_dist 2000 --min_weight 0.5 --out {output} --print_excl ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.excluded_nodes.csv
    """
```

* Make a list of LDpruned SNPs in R

!!! Rule was not executed, script was executed in Rstudio on my local advice (try to implemnt next time)

```
#rule LDpruned_SNPlist:
#  input:
#   'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id'
#   #arg2 = 'angsd/angsd_LC_GL2_cutoff.nInd55.mafs.gz',
#   #arg3 = 'list/LDpruned_snps.list'
#  output:
#    touch('LDpruned_SNPlist.done')
#  log: 'log/LDpruned_SNPlist.log'
#  message:
#    """ Make a list of LDpruned SNPs in R """
#  shell:
#    """
    Rscript scripts/getLDpruned_SNPlist.R
    """
```

### ANGSD: Calculate genotype likelihoods using a list of predefined LD pruned sites

* Index SNP list using angsd

```
rule angsd_index_sites:
  input:
    'list/LDpruned_snps.list'
  output:
    touch('list/index_SNP.done')
  log:
    'log/index_SNP.done'
  threads: 12
  message:
    """ Index SNP list using angsd """
  shell:
    """
    angsd sites index {input} 2> {log}
    """
```

* Get the chromosomes/scaffolds for which we have sites we want to analyse

```

rule getChrom_from_sites:
  input:
    'list/LDpruned_snps.list'
  output: 
    touch('list/getChrom_from_sites.done')
  log:
    'log/getChrom_from_sites.done'
  message:
    """ Get the chromosomes/scaffolds for which we have sites we want to analyse """
  shell:
    """
     cut -f1 {input} | sort | uniq > list/LDpruned_Chrom.list
    """
```

* Calculate genotype likelihoods on predefined LD pruned sites

```

rule angsd_GL2_LDpruned:
  input:
    ref = config["ref"],
    #touched = 'list/index_snp_dgal1.done',
    #touched2 = 'list/getChrom_dgal1_from_sites.done'
    touched = 'list/index_SNP.done',
    touched2 = 'list/getChrom_from_sites.done'
  output:
    touch('angsd/angsd_LC_GL2_maf0018_LDpruned.done')
    #touch('angsd/angsd_LC_GL2_maf0018_LDpruned_dgal1.done')
  log:
    'log/angsd_LC_GL2_maf0018_LDpruned.log'
    #'log/angsd_LC_GL2_maf0018_LDpruned_dgal1.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods on predefined LD pruned sites only using angsd (-doMajorMinor 3 -sites)  """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -anc {input.ref} -out angsd/angsd_LC_GL2_maf0018_LDpruned -GL 2 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doCounts 1 -sites 'list/LDpruned_snps.list' -rf list/LDpruned_Chrom.list 2> {log}
    """
```

#### PCAngsd: Calculate a covariance matrix on pruned data


* estimate covariance matrix

```
rule PCAngsd_LDpruned:
  input:
    touched = 'angsd/angsd_LC_GL2_maf0018_LDpruned.done'
    #GL2 = 'angsd/angsd_LC_GL2_LDpruned.beagle.gz'
  output:
    touch('pcangsd/PCAngsd_GL2_maf0018_LDpruned.done')
  log: 'log/PCAngsd_GL2_maf0018_LDpruned_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_LC_GL2_maf0018_LDpruned.beagle.gz -o pcangsd/PCAngsd_GL2_maf0018_LDpruned_covmat 2> {log}
    """
```

* plot covariance matrix

```
rule plot_covMat_LDpruned:
  input:
    touched = 'pcangsd/PCAngsd_GL2_maf0018_LDpruned.done'
    #covMat = 'pcangsd/PCAngsd_GL2_LDpruned_covmat.cov'
  output:
    pdf = 'pcangsd/PCAngsd_GL2_maf0018_LDpruned.pdf'
  log: 'log/PCAngsd_GL2_LDpruned_plotCovMat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    Rscript scripts/plot_covMat_LDpruned.R pcangsd/PCAngsd_GL2_maf0018_LDpruned_covmat.cov {output.pdf} 2> {log}
    """
```

summary
* angsd GL2, nInd 55, maf 0.018 all sites (angsd_LC_GL2_cutoff_maf0018_nInd55.mafs.gz: 6,430,625 sites)
* angsd GL2, nInd 55, maf 0.018 on LD pruned data (angsd_LC_GL2_maf0018_LDpruned.mafs.gz: 907,194 sites)


Results PCA: Principal-component analysis (PCA) showed three distinct genetic groups matching our reference clonal lines from extant populations D. longispina (blue), D. galeata (orange), D. cucullata (green). PC1 vs PC2 shows that all resting egg samples lie in between D. longispina and D. galeata. 

## Infer admixture proportions using PCAngsd

```
rule PCAngsd_admix:
  input:
    'angsd/angsd_LC_GL2_maf0018_LDpruned.done'
  output:
    touch('pcangsd/PCAngsd_GL2_maf0018_LDpruned_admix_K3.done')
  log: 'log/PCAngsd_GL2_maf0018_LDpruned_admix_K3.log'
  threads: 12 
  message:
    """ Infer admixture proportions using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_LC_GL2_maf0018_LDpruned.beagle.gz -admix -admix_K 3 -o pcangsd/PCAngsd_GL2_maf0018_LDpruned_admix_K3 2> {log}
    """
```

### plot admixture proportions usinfg Rscript

```

setwd("/mnt/data/snakemake_mach2/Daphnia_RestEggs_snakemake_pbs/")

#Load the covariance matrix
admix = npyLoad("pcangsd/PCAngsd_GL2_maf0018_LDpruned_admix_K3.admix.3.Q.npy")
admix

#We will also add a column with population assingments
sample_data <- read.csv("list/metadata_RESTEGGS_LC_OLDvsNEW_badOUT.csv", header = T, sep='\t')
sample_data$sample_id
sample_data$species

# pops
pop <- as.vector(sample_data$species)
pop

# IDs
ID <- as.vector(sample_data$sample_id)
ID



admix.id = as.data.frame(cbind(pop, ID, admix))
names(admix.id) = c("pop","ID","longispina","galeata","cucullata")


# Colour palette with black:
#cbp2 <- c("#009E73", "#000000", "#D55E00", "#56B4E9", "#E69F00", "#F0E442", "#0072B2", "#CC79A7")
# 
pdf("pcangsd/PCangsd_admix_LDpruned_K3_plot.pdf")
par(mar = c(5,4,4,2))
barplot(t(as.matrix(subset(admix.id, select=longispina:cucullata))), col=c("#56B4E9","#D55E00","#009E73"), names=ID, border=NA, las=2, cex.names = 0.5, xlab = "Individuals", ylab = "Ancestry", main = "K3")
dev.off()

jpeg("pcangsd/PCangsd_admix_LDpruned_K3_plot.jpg", height = 431, width = 1057)
par(mar = c(5,4,4,2))
barplot(t(as.matrix(subset(admix.id, select=longispina:cucullata))), col=c("#56B4E9","#D55E00","#009E73"), names=ID, border=NA, las=2, cex.names = 0.5, xlab = "Individuals", ylab = "Ancestry", main = "K3")
dev.off()

png("pop_stats/PCangsd_admix_LDpruned_K3_plot.png", width = 6, height = 6, units = 'in', res = 300)
par(mar = c(5,4,4,6))

barplot(t(as.matrix(subset(admix.id, select=longispina:cucullata))), legend.text = TRUE, args.legend = list(x="right", inset= c(-0.25,0), cex=0.8), col=c("#56B4E9","#D55E00","#009E73"), names=ID, border=NA, las=2, cex.names = 0.5, xlab = "Individuals", ylab = "Ancestry", main = "K3")
dev.off()

```

## Calculate heterozygosity per sample using SFS estimation in angsd

```
rule angsd_saf:
  input:
    ref = config["ref"],
    realigned = 'realigned/{sample}.minq20.realigned.bam'
  output:
    touch('angsd/{sample}.GL2.saf.idx.done')
  log:
    'log/{sample}.GL2.saf.log'
  threads: 24
  message:
    """ SFS Estimation for single samples - part1 """
  shell:
    """
    angsd -i {input.realigned} -out angsd/{wildcards.sample}.GL2 -anc {input.ref} -dosaf 1 -fold 1 -ref {input.ref} -C 50 -minQ 20 -minMapQ 30 -GL 2 2> {log}
    """
```

```
rule real_SFS:
  input:
    'angsd/{sample}.GL2.saf.idx.done'
  output:
    'angsd/{sample}.GL2.est.ml'
  log:
    'log/{sample}.GL2.est.ml.log'
  threads: 24
  message:
    """ SFS Estimation for single samples - part1 """
  shell:
    """
    path=(/home/uibk/c7701178/.conda/envs/eggs/bin/)
    $path/realSFS angsd/{wildcards.sample}.GL2.saf.idx > {output} 2> {log}
    """
```

```
rule plotHeterozygosity:
  input: 'angsd/{sample}.GL2.est.ml'
  output:
    'angsd/{sample}.GL2.heterozygosity.txt'
  log:
    'log/{sample}.GL2.heterozygosity.log'
  message:
    """ plot heterozygosity """
  shell:
    """
   Rscript scripts/plotHeterozygosity.R {input} {output} 2> {log}
    """
```

Rscript to plot heterozygosity


```
#library(tidyverse) #load the tidyverse package for formatting and plotting
library(viridis)
library(plyr)
library(readr)
library(ggplot2)

#In this script we will plot heterozygosity computed by angsd using the -doSaf with the reference genome providing the ancestral state.
#This was achieved by dividing the number of heterozygous sites (the second entry in the SFS)
#by the total number of sites (first entry + second entry in the SFS)
#to obtain the proportion of heterozygous sites for each genome.


setwd("/mnt/data/snakemake_mach2/Daphnia_RestEggs_snakemake_pbs/")

#args <- commandArgs(trailingOnly=TRUE)

# population assingments
sample_data <- read.csv("list/metadata_RESTEGGS_LC_OLDvsNEW_badOUT.csv", header = T, sep='\t')
sample_data$sample_id
sample_data$species

# species
species <- as.vector(sample_data$species)
species

# IDs
ID <- as.vector(sample_data$sample_id)
ID


basedir <- "angsd" # Make sure to edit this to match your $BASEDIR
bam_list <- list.files(path = basedir, pattern = "*est.ml", full.names = TRUE)

bam_list



data_tsv = ldply(bam_list, scan)

# total number of sites
totalNR_sites <- (data_tsv$V1 + data_tsv$V2)
totalNR_sites

# number of homozygote sites
NR_HOMsites <- data_tsv[1]

# number of heterozygote sites
NR_HETsites <- data_tsv[2]
NR_HETsites

# individual heterozygosity
IND_heterozygosity <- (NR_HETsites/totalNR_sites)

output <- cbind(ID, species, NR_HOMsites, NR_HETsites, IND_heterozygosity)
output

names(output) = c("ID", "species", "NR_HOMsites", "NR_HETsites", "IND_heterozygosity")
options(scipen = 100)
write.table(output,"pop_stats/heterozygosity_angsd.txt", sep ="\t", quote = F)


# plot bar charts of heterozygosity estimates

barplot(t(as.matrix(output$IND_heterozygosity)), names=output$ID, border=NA, las=2, cex.names = 0.5, cex.axis = 0.8, xlab = "Individuals", ylab = "Heterozygosity")


png("pop_stats/heterozygosity_angsd.png", width = 4, height = 4, units = 'in', res = 300)
# barplot
p<-ggplot(data=output, aes(x=ID, y=IND_heterozygosity, color=species, fill=species, label=ID))+
  geom_bar(stat="identity", width=0.5)+
    theme(axis.text = element_text(size = 6))+
      scale_color_manual(values=c("#009E73", "#000000", "#D55E00", "#56B4E9"))
      p2 <- p + scale_fill_manual(values=c("#009E73", "#000000", "#D55E00", "#56B4E9"))+
        xlab("Indiviuals")+
	  ylab("Heterozygosity")+
	    coord_flip()
	    p2 + theme(panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
	      scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
	      dev.off()


```


