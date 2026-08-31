## Lab # 7 RNA-Seq

## Table of Contents
1. [Introduction](#intro)
2. [Setting Up](#setup)
3. [The Experiment](#experiment)
4. [Data Cleaning](#clean)
5. [Mapping Reads to the Human Genome](#mapping)
6. [Transcript Assembly](#transcripts)
7. [Differential Gene Expression Analysis](#dge)
8. [Interpretation](#interpretation)

<a name="intro"></a>
## Introduction

The goal of the lab is to introduce bulk RNA-sequencing analysis on the command line. The lab will walk through quality control and alignment on the command line along with analysis and data visualization implemented in R.

**Flash Updates**
* *RNA-Seq* 
* *False Discovery Rate (FDR)* 
* *Principal Components Analysis (PCA)* 

**Demo Videos**
* [Logging onto the Cluster](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/houlahke_mcmaster_ca/IQC2ftW_ZyPVS5qzMttT4FYvAX0pB_x6evXgf2O9YmDjYWI?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=WCQxQd) ~2 minutes
* [Viewing FASTQC Results](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/houlahke_mcmaster_ca/IQCDQp8AuxWwTbQblBHpQmidAcPwmwdYdprqX7dhEixl2Fk?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=fNC7WY) ~3 minutes

**Computer Resources**
* The lab can be completed on the *cluster*. Refer back to Lab 1 and 2 to refresh using command line arguments in Linux and R.

**Grading**
* Questions are for your learning and are not graded
* Problems are worth 5 points each (-1 for each error)
* Submit your answers to the Problems, plus any supplmental multiple choice questions, on **A2L Quizzes** before the deadline
* An answer key to Questions and Problems will be provided on A2L after the deadline

<a name="qc"></a>
## Setting Up

Today’s lab will use the cluster. **Reminder:** to log onto the cluster you must be on the McMaster network or VPN. 

```bash
ssh -l <macid> acf-access-student.csu.mcmaster.ca
```

**Annotation (GTF format) File(s):**

All of the below files can be found `/workspace/lab/studentlab/lab7_rnaseq/references`

```Bash
gencode.v29.annotation.gtf.gz
```

**Illumina Sequencing (FASTQ format) Files(s):**

All of the below files can be found `/workspace/lab/studentlab/lab7_rnaseq/data`

```
adrenal.fastq
HLE_Cd_1_forward.fastq.gz
HLE_Cd_1_reverse.fastq.gz
HLE_Cd_2_forward.fastq.gz
HLE_Cd_2_reverse.fastq.gz
HLE_Cd_3_forward.fastq.gz
HLE_Cd_3_reverse.fastq.gz
HLE_Ctrl_1_forward.fastq.gz
HLE_Ctrl_1_reverse.fastq.gz
HLE_Ctrl_2_forward.fastq.gz
HLE_Ctrl_2_reverse.fastq.gz
HLE_Ctrl_3_forward.fastq.gz
HLE_Ctrl_3_reverse.fastq.gz
```

<a name="experiment"></a>
## The Experiment

> Flash Update - RNA-Seq

We are going to examine the response of the human transcriptome in a human lens epithelial cell line (part of the eye) exposed to Cadmium, as preliminary microarray work has suggested Cadmium exposure, via the MTF-1 transcription factor, impacts lens development and maintenance. The experiment is RNA-Seq of three Cadmium exposed replicates and 3 Control replicates, using the GRCh38 version of the human genome annotation as reference. The RNA-Seq was performed using an Illumina HiSeq with 2 x 50 bp mate pair sequencing.

We are going to start by using the FastQC tool to examine the quality of some of the RNA-Seq data. Details on all the plots can be found here: [FASTQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) and [video tutorial](http://www.youtube.com/watch?v=bz93ReOv87Y).

```bash
# create directory for lab 7
mkdir lab7_rnaseq
cd lab7_rnaseq
# let's make symlinks to the file in our directory
ln -s /workspace/lab/studentlab/lab7_rnaseq/data/HLE_Cd_1_forward.fastq.gz .
# run fastqc on foward fastq
fastqc HLE_Cd_1_forward.fastq.gz
# create directory to store initial fastqc run
mkdir raw_fastqc
# move output into folder
mv HLE_Cd_1_forward_fastqc.html raw_fastqc/
mv HLE_Cd_1_forward_fastqc.zip raw_fastqc/
```

**Question #1. How many mRNA were sequenced from each replicate and does this data need any adaptor removal or quality trimming?**

**Problem #1. This lab is using only a fraction of the total data so it does not take too long, but also perform FASTQC on the full replicate from a different experiment (adrenal.fastq). When a full RNA-Seq run is analyzed, do the samples pass FastQC's quality control checks for per-sequence GC content and sequence duplication levels? These checks passed when we were assembling bacterial genomes. If they do not pass for these RNA-Seq data, suggest reasons.**

<a name="clean"></a>
## Data Cleaning

Before we start to compute on the data, let's make sure we are on a compute node.

```Bash
# request an interactive node on the cluster with 2G that will run for 24 hours
salloc --mem=2G --time=24:00:00 --partition=classroom
srun --pty bash
```

Even if the data as a whole passed FASTQC, quality trimming and filtering is still highly recommended to remove or trim individual sequences of poor quality. Run *Trimmomatic* (paired-end with separate input files, plus ILLUMINACLIP with TruSeq3 for paired-end MiSeq or HiSeq) on all the samples. For example:

```bash
# create directory to store trimmed data in
mkdir trimmomatic
# run trimmomatic in paired-end mode
apptainer run -B /workspace/lab:/workspace/lab /workspace/lab/studentlab/lab5_genome_assembly/trimmomatic_v0.40.sif trimmomatic PE -threads 1 \
    HLE_Cd_1_forward.fastq.gz \
    HLE_Cd_1_reverse.fastq.gz \
    trimmomatic/HLE_Cd_1_forward_trimmed_paired.fq.gz \ 
    trimmomatic/HLE_Cd_1_forward_trimmed_unpaired.fq.gz \ 
    trimmomatic/HLE_Cd_1_reverse_trimmed_paired.fq.gz \ 
    trimmomatic/HLE_Cd_1_reverse_trimmed_unpaired.fq.gz \ 
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:True \
    SLIDINGWINDOW:4:20 
```

`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10` removes the adapters for Illumina TruSeq3 sequencing which is the method used to generate these data. 

`SLIDINGWINDOW:4:15` scans the read with a 4-base wide sliding window, cutting where the average quality per base drops below 15.

**Note**: Trimmomatic under these settings creates both **paired** and **unpaired** output. We only want to use paired reads in our data, so will ignore the unpaired files.

**Question #2. Run *FASTQC* on a couple of your samples to see if the data has changed in quality. Has anything improved?**

<a name="mapping"></a>
## Mapping Reads to the Human Genome

Before we can interpret these data, we need to map the FASTQ reads to the reference human genome (hg38). We cannot use the standard Burrows-Wheeler Transform software BWA or Bowtie, since RNA-Seq data needs to be corrected for introns and exons. Instead, we will use [HISAT2](https://daehwankimlab.github.io/hisat2/main/), which can handle splice junction boundaries as well as control for fragment sizes.

The first step is download the [reference files](https://daehwankimlab.github.io/hisat2/download/) required by *HISAT2*. This step has already been done for you. The data can be found `/workspace/lab/studentlab/lab7_rnaseq/references/grch38`.

Perform *HISAT2* read mapping for each sample, using the reference files you just downloaded. The command is provided below, however, please submit as a job to slurm. **Do not run on head node.** It is recommended to run the alignment with 6GB of memory requested. As a reminder of how to create a *sbatch script* please refer back to Lab 1. Also specifically for hisat2 you will have to load the module in advance using the command `module load hisat2_2.2.3` in your sbatch script.

```bash
mkdir hisat2
# for hisat2 we need to load it 
module load hisat2_2.2.3
# run HISAT2 to align reads to reference genome
hisat2 --add-chrname -x /workspace/lab/studentlab/lab7_rnaseq/references/grch38/genome \
-1 trimmomatic/HLE_Cd_1_forward_trimmed_paired.fq.gz \
-2 trimmomatic/HLE_Cd_1_reverse_trimmed_paired.fq.gz \
-S hisat2/HLE_Cd_1.sam
```

The [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) output is standard file format for sequencing alignment. However, the SAM file format is uncompressed and can take up a lot of space. Thus, we can convert SAM files to BAM files which are the binary equivalent.

```bash
# convert sam file to bam file
samtools view -bS hisat2/HLE_Cd_1.sam > hisat2/HLE_Cd_1.bam
# once the file is coverted, we can remove the original sam file
rm hisat2/HLE_Cd_1.sam
```

Finally, we need to assess our alignment. By default, *HISAT2* provides information on alignment to the *standard out*.

**Problem #2. HISAT2 creates a BAM file that contains the alignment information. What percentage of read pairs aligned uniquely to one location in the genome and what percentage may represent multiple copy genes? What was the overall alignment rate? Would you say this is a good RNA-Seq data set? Why?**

<a name="transcripts"></a>
## Transcript Assembly

Now that the raw RNA-Seq data have been aligned to the reference human genome, we can assemble the data into individual transcripts as a step towards identifying differential gene expression (DGE). The *htseq-count* tool determines the transcripts at each gene in the reference and provides un-normalized counts.

Perform *featureCounts* on each replicate's *HISAT2* BAM file, using the *gencode.v29.annotation.gtf.gz* annotation file and the *Reverse* stranded option (which reflects use of a first-strand synthesis kit during library construction, see [PMID 32415774](https://pubmed.ncbi.nlm.nih.gov/32415774/)). Reverse strandedness can be implemented in featureCounts using the parameter `-s 2`. For example:

```bash
# create new output directory
mkdir feature_counts
# run feature counts 
apptainer run -B /workspace/lab:/workspace/lab /workspace/lab/studentlab/lab7_rnaseq/scripts/feature-counts_v2.0.0.sif featureCounts -p -a /workspace/lab/studentlab/lab7_rnaseq/references/gencode.v29.annotation.gtf.gz \
-o feature_counts/HLA_Cd_1_feature_counts.txt \
-s 2 \
hisat2/HLE_Cd_1.bam
```

**Question #3. Examine the results of featureCounts. How many total transcripts are quantified? Write a simple R script to determine how many transcript have 3 or more reads mapped to them. What proportion of transcripts have at least 3 reads?**

Repeat the above steps for the remaining samples.

<a name="dge"></a>
## Differential Gene Expression Analysis

> Flash Update - False Discovery Rate (FDR)

We are going to use *DESeq2* to both normalize and perform significance tests on these data. To do this, we can run the script `run_deseq2.R`.

```bash
# run DESeq2
Rscript /workspace/lab/studentlab/lab7_rnaseq/scripts/run_deseq2.R
```

*DESeq2* will create a results file that included significance testing (using the P-adj to reflect correction for false discovery), a principal components plot to visualize differences in overall transcriptome among the replicates, and a table of normalized counts.

> Flash Update - Principal Component Analysis (PCA)

**Question #4. Look at transcript differential expression testing and then try filter in R for significant differences in transcript abundance (P-adj < 0.05). How many genes are differentially expressed in this experiment at this corrected alpha value?**

**Question #5. Look at the normalized counts and then try sort in R to determine the most highly expressed gene in Cadmium exposed cells. Is it the same for each replicate?**

<a name="interpretation"></a>
## Interpretation

At this point, we have a robust statistical analysis of these RNA-Seq data, with a resulting list of significantly differentially expressed genes, that are labeled using *ENSEMBL_GENE_ID* identifiers. We will be using [gProfiler](https://biit.cs.ut.ee/gprofiler/gost) to provide some biological context. gProfiler identifies biological pathways that are enriched more than expected amongst the list of differentially expressed genes.

Take your list of differentially abundance transcript (Question 4) and paste them into the input box on the gProfiler page. Set the parameters according to the screenshot below. In this case, we will just be considering the KEGG pathways and pathways will be considered significant if FDR < 0.05.

![gProfiler](gprofiler.png)

**Problem #3. What is your overall interpretation of the impact of Cadmium on human lens epithelial cells?**

