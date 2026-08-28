## Lab # 6 - Genome alignment and variant detection

## Table of Contents
1. [Introduction](#intro)
2. [Data](#data)
3. [Quality control](#quality)
4. [Alignment](#alignment)
5. [Variant Calling](#variant)
6. [Annotation](#annotation)

<a name="intro"></a>
## Introduction

The goal of this lab is to walk through a typical DNA alignment and variant detection workflow.

**Flash Updates**
* *single nucleotide polymorphisms (SNPs)* 
* *sequencing alignment map (SAM/BAM)*
* *variant calling format (VCF)*

**Demo Videos**
* [Logging onto the Cluster](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/houlahke_mcmaster_ca/IQC2ftW_ZyPVS5qzMttT4FYvAX0pB_x6evXgf2O9YmDjYWI?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=WCQxQd) ~2 minutes
* [Viewing FASTQC Results](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/houlahke_mcmaster_ca/IQCDQp8AuxWwTbQblBHpQmidAcPwmwdYdprqX7dhEixl2Fk?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=fNC7WY) ~3 minutes
* [Submitting to Slurm](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/houlahke_mcmaster_ca/IQAmpCiUipYAQ4oNEQRFXhLxAe2xC2hnFuVGnkMGDs1Y7eM?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=diY7N6) ~5 minutes

**Computer Resources**
* Today’s lab will use the cluster. **Reminder:** to log onto the cluster you must be on the McMaster network or VPN. 

```bash
ssh -l <macid> acf-access-student.csu.mcmaster.ca
```

**Grading**
* Questions are for your learning and are not graded
* Problems are worth 5 points each (-1 for each error)
* Submit your answers to the Problems, plus any supplmental multiple choice questions, on **A2L Quizzes** before the deadline
* An answer key to Questions and Problems will be provided on A2L after the deadline

<a name="data"></a>
## Data

> Flash Update - SNPs

In the previous lab you assembled a *Salmonella enterica* genome sequence, learning the steps of how to filter and assemble raw Illumina sequencing reads to form genome contigs and scaffolds. Now you are going to explore a human genome sequence and, rather than assemble it *de novo*, you are going to align the DNA sequences to a reference genome and detect SNPs.

We are going to be working with whole exome sequencing (WXS) of a breast cancer cell line (MDA-MB-415). Because we are only sequencing the exome, this data represents 2% of the full human genome. 

We will be using the below data: 
```Bash
/workspace/lab/studentlab/lab6_genome_variants/SRR8619134_1.fastq.gz
/workspace/lab/studentlab/lab6_genome_variants/SRR8619134_2.fastq.gz
```

Again we can create a new directory in our home directory and create simlinks to these FASTQ files.

```Bash
# create new directory for lab 5
cd ~
mkdir lab6_genome_variants
cd lab6_genome_variants
# create simlinks
ln -s /workspace/lab/studentlab/lab6_genome_variants/SRR8619134_1.fastq.gz SRR8619134_1.fastq.gz
ls -s /workspace/lab/studentlab/lab6_genome_variants/SRR8619134_2.fastq.gz SRR8619134_2.fastq.gz
```

As before, run the below commands off the symlinks.

<a name="quality"></a>
## Quality control

Before we start to compute on the data, let's make sure we are on a compute node. We can do this by logging onto an interactive node as we did in lab 2 and 5. 


```Bash
# request an interactive node on the cluster with 2G that will run for 24 hours
salloc --mem=2G --time=24:00:00 --partition=classroom
srun --pty bash
```

As we did in lab 5, we are going to use FASTQC to determine the quality of our data. Details on all the plots can be found here: [FASTQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) and [video tutorial](http://www.youtube.com/watch?v=bz93ReOv87Y). See the demo video for how to view the FASTQC results.

```bash
# run fastqc on original fastq files
fastqc SRR8619134_1.fastq.gz
fastqc SRR8619134_2.fastq.gz
```

**Question #1. What is the quality of your data? Is there any evidence that the sequence library is biased (i.e. non-random)? Explain your reasoning.**

**Question #2. Is there any evidence that the reads still have their adapters? Do we need to trim the adapters?**

<a name="alignment"></a>
## Alignment 

> Flash Update - SAM/BAM

Next, we need to figure out where each sequencing read comes from in the genome. To do this, we are going to align our sequencing reads to the human reference genome. In this case, we are going to use [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/).

To align the reads to the reference, we are going to use `BWA MEM`. The reference and its indices can be found `/workspace/lab/studentlab/lab6_genome_variants/GRCh38/GRCh38.d1.vd1.fa`. This is going to take a while so it is highly recommended to submit this command to Slurm by creating a sbatch script (see demo video).

```Bash
# run alignment using bwa mem
bwa mem /workspace/lab/studentlab/lab6_genome_variants/GRCh38/GRCh38.d1.vd1.fa SRR8619134_1.fastq.gz SRR8619134_2.fastq.gz > SRR8619134.sam
```

By default, `BWA MEM` takes the fastq files and outputs a `SAM` (<u>S</u>equence <u>A</u>lignment <u>M</u>ap) file. More details on the format of a SAM file can be found [here](https://samtools.github.io/hts-specs/SAMv1.pdf). SAM files are large and take up a lot of space. Therefore, to save space, we can convert SAM files into their binary equivalent, called a `BAM` file.

```Bash
# sort a SAM file and convert it to a BAM file 
samtools sort SRR8619134.sam -O BAM -o SRR8619134.sort.bam
# index bam file 
samtools index SRR8619134.sort.bam
```

samtools fastq -1 output_R1.fastq.gz -2 output_R2.fastq.gz \
    -0 /dev/null -s /dev/null -n sorted_by_name.bam

**Question #3. What is the file size of the original SAM file? How does it compare to the binary BAM file?**

Once you have your `BAM` you can delete your `SAM` file to save storage space. They both include the same data but the `BAM` file is more compressed.

```Bash
# remove SAM file
rm SRR8619134.sam
```

Let's look to see how well our data aligned to the reference genome using `flagstat`. To understand how to interpret the output the flagstat command, check out the [documentation](https://www.htslib.org/doc/samtools-flagstat.html).

```Bash
samtools flagstat SRR8619134.sort.bam
```

**Problem #1. Do you have a high quality alignment? Consider what proportion of reads are mapping and whether they are mapping in the expected alignment with their mate.**

**Problem #2. Why do not all reads align? Provide some possible explanations.**

<a name="variant"></a>
## Variant Calling

> Flash Update - VCF 

Next, we want to find DNA variants. Here we are going to focus on SNPs. We are going to use `bcftools`. The first step is to run `mpileup` which calculates the number of reads and which base are found at each site. Then, `bcftools call` uses this information to find positions where the base detected differs from the reference genome.

```Bash
# run bcftools mpileup and call to detect possible SNPs
bcftools mpileup \
    -f /workspace/lab/studentlab/lab6_genome_variants/GRCh38/GRCh38.d1.vd1.fa \
    SRR8619134.sort.bam \
    | bcftools call -mv -Oz \
    -o variants.vcf.gz
```

```Bash
bcftools stats variants.vcf.gz > variants_stats.txt
```

**Question #4. How many SNPs did you detect? Is this what you would expect? What evidence supports each variant?**

Not all variants may be real SNPS. We often need to apply additional filtering to reduce false positives in our calls. There are a few criteria we can implement to enrich for high confidence SNP calls:
* Filter out SNPs with insufficient sequencing reads at the site (e.g. `INFO/DP<10`)
* Filter out SNPs at sites with low mapping quality at site (e.g. `INFO/MQ<40`)
* Filter out SNPs with low genotyping quality (e.g. `QUAL<30`)

```Bash
# index variants
tabix variants.vcf.gz 
# filter variants
bcftools view \
    -e 'INFO/DP<10 || INFO/MQ<40 || QUAL<30' \ 
    -v snps \
    -R /workspace/lab/studentlab/lab6_genome_variants/exome.bed \
    -o variants_filtered.vcf.gz \
    -O z \
    variants.vcf.gz
```

**Question #5. Why is it important to include the `-R exome.bed` parameter?**

**Question #6. How many SNPs are detected after filtering? How many were filtered out by our depth filter? By our quality filter? By the genotyping confidence filter?**

<a name="annotation"></a>
## Annotation

There are multiple tools we can use to annotate our SNPs to determine which genes they may be falling in and what their predicted function may. 

```Bash
# annotate variants
java -jar /opt/COMMON_APPLICATIONS/snpEff/snpEff.jar GRCh38.mane.1.0.ensembl variants_filtered.vcf.gz > variants_annotated.vcf
```

We can reformat the vcf to extract gene names and extract all variants predicted to be high impact.

```Bash
# extract high impact variants 
bcftools view -i 'INFO/ANN ~ "HIGH"' -O z -o high_impact.vcf.gz variants_annotated.vcf
# reformat output
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
    high_impact.vcf.gz | \
    awk -F'\t' '{split($5,a,"|"); print $1,$2,$3,$4,a[4]}' OFS='\t' > high_impact_genes.txt
```

**Problem #3. What variant(s) could be driving the tumour that the cell line was derived from? Justify your answer.**

You may also want to consider `/workspace/lab/studentlab/lab6_genome_variants/Cosmic_breast_genes.txt` which provides a list of possible driver genes in breast cancer. 