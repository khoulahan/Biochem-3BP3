## Lab # 8 DNA Methylation

## Table of Contents
1. [Introduction](#intro)
2. [Setting Up](#setup)
3. [The Experiment](#experiment)
4. [Downloading data from GEO](#data)
5. [Cleaning the data](#clean)
6. [Data visualization](#visualization)
7. [Identifying differential methylation](#stats)

<a name="intro"></a>
## Introduction

The goal of this lab is to introduce you to DNA methylation data. This lab will have you download public data from the Gene Expression Omnibus (GEO), understand the data type you download, run some basic linear models, and explore different data visualization in ggplot2. 

We will be looking at 450K microarray data, as it's less computationally intensive. This tutorial should be able to be completed on your laptop, however, the instructions below walk through how to run it on the cluster.

**Flash Updates**
* *CpG Islands* 
* *Bisulfite Conversion* 
* *DMRs* 

**Demo Videos**
* [Logging onto the Cluster](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/houlahke_mcmaster_ca/IQC2ftW_ZyPVS5qzMttT4FYvAX0pB_x6evXgf2O9YmDjYWI?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=WCQxQd) ~2 minutes

**Computer Resources**
* The lab can be completed on the *cluster*. Refer back to Lab 1 and 2 to refresh using command line arguments in Linux and R.

**Lab Resources**

*Relevant publications:* GSE100197
https://academic.oup.com/hmg/article/27/1/135/4581994?login=true

*Relevant educational videos:*
*DreamRs, an introductory interface to ggplot2:*
https://dreamrs.shinyapps.io/esquisse/

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

Before we start to compute on the data, let's make sure we are on a compute node. For this lab we are going to need to request 8G.

```Bash
# request an interactive node on the cluster with 10G that will run for 24 hours
salloc --mem=10G --time=24:00:00 --partition=classroom
```

```Bash
# log onto the interactive node
srun --pty bash
```

<a name="experiment"></a>
## The Experiment

> Flash Update - CpG Islands

In this lab we will assess placental DNA methylation changes in early- and late- onset preeclampsia, a maternal hypertensive disorder during pregnancy. This lab will walk you through downloading the DNA methylation data and associated metadata, cleaning the data, a linear regression, and example data visualization with ggplot2. Following this lab, I want each of you to take the data and be creative in using ggplot2 (either in R or dreamRs interface) to show something different with this data and explain how you would interpret your results.

<a name="download"></a>
## Microarray data: Downloading data from GEO

Let's create a new folder to store all our outputs for lab 8.

```bash
# make sure you are in your home directory
cd ~
# make lab 9 directory
mkdir lab8_methylation
# enter new directory
cd lab8_methylation
```

Now we are going to download DNA methylation data from a database known as Gene Expression Omnibus (GEO). Every dataset on GEO has a corresponding accession ID used to identify it. So we are going to download the dataset corresponding to accession ID `GSE100197`. To do so, run the following:

```Bash
# download DNA methylation data and metadata 
Rscript /workspace/lab/studentlab/lab8_methylation/download_dna_methylation_data.R
```
You should now have the following files in your home directory: 
1. `GSE100197_methylation_data.txt` (reformatted methylation data)
2. `GSE100197_clinical_metadata.csv` (reformatted metadata)
3. `GSE100197_series_matrix.txt.g`z (raw methlyation data downloaded from GEO)
4. `GSE100197_family.soft.gz` (raw metadata downloaded from GEO)

We are going to read in the methylation data and inspect it. Load an R session: 

```R
# load R packages
library(rlang)
library(limma)
library(ggplot2)

# read in dna methylation data 
data <- read.delim('GSE100197_methylation_data.txt', as.is = TRUE)

# Inspect
class(data)
dim(data) # What does this tell you about the data?
head(data)
str(data) # What does this tell you about the data?

class(data$ID_REF)
data$ID_REF <- as.factor(data$ID_REF) #Why am I doing this?
```

**Question #1. What does `dim(data)` tell you about the data?**

**Question #2. What does `str(data)` tell you about the data?**

**Question #3. How do you change how `ID_REF` is coded? Why do you need to change how `ID_REF` is coded?**

**Problem #1. What is the difference between `chr`, `factor`, and `num`? Why do variables need to be set to a specific variable type?**

Next, we are going to inspect the metadata. 

```R
# read in metadata 
metadata <- read.csv('GSE100197_clinical_metadata.csv')

# Inspect available clinical variables
print(dim(metadata))
print(names(metadata))
print(head(metadata[, seq_len(min(8, ncol(metadata)))]))
print(dim(metadata))
print(str(metadata))
```
**Question #4. Compare the information you get from looking at data details with information you get from looking at metadata details.**

**Question #5. What does `dim(metadata)` tell you about the data?**

**Question #6. What does `str(metadata)` tell you about the data? What do the columns and data mean?**


<a name="clean"></a>
## Microarray data: Cleaning the data

> Flash Update - Bisulfite conversion

Next we are going to clean up the metadata structure and check that all samples match between the methylation data and the metadata

```R
#Let's fix the metadata structure
# We only need a subset of the metadata
metadat <- metadata[,c(1,5:12)]
str(metadat)
metadat$geo_accession <- as.factor(metadat$geo_accession)
metadat$subject.id <- as.factor(metadat$subject.id)
metadat$pathology.group <- as.factor(metadat$pathology.group)
metadat$fetal.sex <- as.factor(metadat$fetal.sex)
metadat$gestational.age <- as.numeric(metadat$gestational.age)
metadat$sample.tissue <- as.factor(metadat$sample.tissue)
metadat$X450k.plate <- as.factor(metadat$X450k.plate) #Denotes the 96 well plate samples were prepped on
metadat$X450k.sentrix_id <- as.factor(metadat$X450k.sentrix_id) # Denotes the microarray chip samples were placed on
metadat$X450k.sentrix_position <- as.factor(metadat$X450k.sentrix_position) #Denotes the row and column samples were placed on the chip

#Let's change data ID_REF to the actual sample names
rownames(data) <- data$ID_REF
data$ID_REF <- NULL
head(data)

#Check to make sure all samples match between data and metadata. This function should return 0 if there are no samples that differ between the data and metadata
setdiff(colnames(data), metadat$geo_accession)
```

Next we are going to identify any replicate samples in both methylation data and teh metadata and remove them.

```R
# Remove replicate samples in both data and metadata
# Identify replicate samples (ignoring case and surrounding whitespace)
is_replicate <- grepl(
  "REPLICATE",
  trimws(as.character(metadata$pathology.group)),
  ignore.case = TRUE
)

# Check which samples will be removed
print(metadata[is_replicate, c("geo_accession", "pathology.group")])

# Get their GEO accessions before removing them
replicate_samples <- metadata$geo_accession[is_replicate]

# Remove replicate samples from metadata
metadat <- metadat[!is_replicate, , drop = FALSE]

# Remove matching sample columns from data
data <- data[, !colnames(data) %in% replicate_samples, drop = FALSE]

#checking to make sure there are the same number of samples in metadata and data
nrow(metadat) == ncol(data)
```

<a name="visualization"></a>
## Microarray data: Data visualization

We are going to use the R package `ggplot2` to generate a PCA plot. PCA stands for principal component analysis which is a dimensionality reduction technique that identifies the axes that can explain the largest amount of variability in your data. PCA was a Flash Update topic in lab 7. 

```R
# Transpose data so rows = samples and columns = features
data_t <- t(as.matrix(data))

# Perform PCA, scaling features to comparable variance
pca <- prcomp(data_t, center = TRUE, scale. = TRUE)

# Create plotting data frame
plot_data <- data.frame(
  geo_accession = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# Match sample metadata using GEO accession
plot_data$sentrix_position <- metadat$X450k.sentrix_position[
  match(plot_data$geo_accession, metadat$geo_accession)
]
```

We are going to make a series of PCA plots. While the x and y-axes will remain the same in each plot, the colouring schemes of the points will change. In the PCA plots, each point represents a sample. PC1 (x-axis) and PC2 (y-axis) are the two axes the explain the most variability in the methylation patterns across the samples. If methylation patterns are similar, the samples should cluster closer to together on the PCA plot.

For each plot, the data points will be coloured by different features of the data to determine how much each feature can explain the variability we observe in the methylation patterns. For starters, we are going to look at the **Sentrix position**. The Sentrix position tells us where on the chip the sample has hybridized. This will allow us to determine if there are any biases by location on chip.

```R
# Plot PCA clustering: Sentrix position
png(filename = "2026_pca_sentrixpos.png", height = 6, width = 6, units = "in", res = 300)
ggplot(plot_data, aes(x = PC1, y = PC2,
                      colour = factor(sentrix_position))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(
    aes(label = sentrix_position),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  labs(
    title = "PCA of samples",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)"),
    colour = "Sentrix position"
  ) +
  theme_bw()
dev.off()
```

If you are running this lab on the cluster, in order to look at the plot you just created you will need to copy it to your local machine. To do so, in the terminal on your local machine run: 

```bash
scp macid@acf-access-student.csu.mcmaster.ca:/workspace/lab/macid/lab8_methylation/2026_pca_sentrixpos.png .
```

Next we are going to look at the **Sentrix ID**. The Sentrix ID tells us which chip the sample has hybridized. This will allow us to determine if there are any biases by chip.

```R
# Plot PCA clustering: SentrixID
plot_data <- data.frame(
  geo_accession = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

plot_data$sentrix_id <- metadat$X450k.sentrix_id[
  match(plot_data$geo_accession, metadat$geo_accession)
]

png(filename = "2026_pca_sentrixID.png", height = 6, width = 6, units = "in", res = 300)
ggplot(plot_data, aes(x = PC1, y = PC2,
                      colour = factor(sentrix_id))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(
    aes(label = sentrix_id),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  labs(
    title = "PCA of samples",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)"),
    colour = "sentrix_id"
  ) +
  theme_bw()
dev.off()
```

Third, we are going to investigate variability by **fetal sex** to determine the methylation profile differs by sex. 

```R
# Plot PCA clustering: Fetal sex
png(filename = "2026_pca_fetalsex.png", height = 6, width = 6, units = "in", res = 300)
plot_data <- data.frame(
  geo_accession = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

plot_data$fetal_sex <- metadat$fetal.sex[
  match(plot_data$geo_accession, metadat$geo_accession)
]

ggplot(plot_data, aes(x = PC1, y = PC2,
                      colour = factor(fetal_sex))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(
    aes(label = fetal_sex),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  labs(
    title = "PCA of samples",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)"),
    colour = "Fetal Sex"
  ) +
  theme_bw()
dev.off()
```

Fourth, we can investigate if methlyation profiles differ by pathology group. 

```R
# Plot PCA clustering: Pathology Group
plot_data <- data.frame(
  geo_accession = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

plot_data$path <- metadat$pathology.group[
  match(plot_data$geo_accession, metadat$geo_accession)
]

png(filename = "2026_pca_pathology.png", height = 6, width = 6, units = "in", res = 300)
ggplot(plot_data, aes(x = PC1, y = PC2,
                      colour = factor(path))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(
    aes(label = path),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  labs(
    title = "PCA of samples",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)"),
    colour = "Pathology"
  ) +
  theme_bw()
dev.off()
```

Finally, we are going to investigate if methylation profiles differ by **gestational age**. 

```R
# Plot PCA clustering: Gestational age
plot_data <- data.frame(
  geo_accession = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

plot_data$gestational_age <- metadat$gestational.age[
  match(plot_data$geo_accession, metadat$geo_accession)
]

plot_data$gestational_age <- as.numeric(as.character(plot_data$gestational_age))

png(filename = "2026_pca_GA.png", height = 6, width = 6, units = "in", res = 300)
ggplot(plot_data, aes(x = PC1, y = PC2, colour = gestational_age)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(
    aes(label = gestational_age),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  scale_colour_gradient(
    low = "blue",
    high = "red",
    name = "Gestational age"
  ) +
  labs(
    title = "PCA of samples",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)")
  ) +
  theme_bw()
dev.off()
```

**Problem #2. What do these plots tell you about how the methylation profiles of these samples differ by Sentrix position, Sentrix ID, fetal sex, pathology group and gestational age?**

**Question #7. What variables do you think should be included in a model to identify differential DNA methylation?**

<a name="stats"></a>
## Microarray data: Identifying differential methylation

> Flash Update - Differentially methylated regions (DMRs)

Finally, we are going to identify differentially methylation regions using an approach similar to what we did in lab 7 (RNA-sequencing). To do so, we are going to use an R package called `limma`.

```R
# Define the pathology groups to include
groups <- c("EOPE", "PreT", "LOPE", "Term")

# Match metadata to the sample columns in data
sample_cols <- intersect(colnames(data), metadat$geo_accession)

meta <- metadat[
  match(sample_cols, metadat$geo_accession),
  ,
  drop = FALSE
]

# Keep only the four pathology groups
keep <- meta$pathology.group %in% groups

meta <- meta[keep, , drop = FALSE]
sample_cols <- sample_cols[keep]

# Extract the numeric data matrix
mat <- as.matrix(data[, sample_cols, drop = FALSE])
storage.mode(mat) <- "numeric"

# Set up covariates
meta$pathology.group <- factor(
  meta$pathology.group,
  levels = groups
)

meta$gestational.age <- as.numeric(
  as.character(meta$gestational.age)
)

meta$fetal.sex <- factor(meta$fetal.sex)
meta$X450k.sentrix_id <- factor(meta$X450k.sentrix_id)

# Confirm group sample sizes
print(table(meta$pathology.group))

# Construct design matrix
design <- model.matrix(
  ~ 0 + pathology.group +
    gestational.age +
    fetal.sex +
    X450k.sentrix_id,
  data = meta
)

colnames(design)

# Fit linear model for every probe
fit <- lmFit(mat, design)

# Define pairwise contrasts
contrast_matrix <- makeContrasts(
  EOPE_vs_PreT = pathology.groupEOPE - pathology.groupPreT,
  EOPE_vs_LOPE = pathology.groupEOPE - pathology.groupLOPE,
  EOPE_vs_Term = pathology.groupEOPE - pathology.groupTerm,
  PreT_vs_LOPE = pathology.groupPreT - pathology.groupLOPE,
  PreT_vs_Term = pathology.groupPreT - pathology.groupTerm,
  LOPE_vs_Term = pathology.groupLOPE - pathology.groupTerm,
  levels = design
)

# Apply contrasts and empirical Bayes moderation
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_contrasts <- eBayes(fit_contrasts)

# Extract results for each comparison
results <- lapply(
  colnames(contrast_matrix),
  function(contrast) {
    topTable(
      fit_contrasts,
      coef = contrast,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
  }
)

names(results) <- colnames(contrast_matrix)

# Example: inspect EOPE versus Term
head(results$EOPE_vs_Term)

# Export results
for (contrast in names(results)) {
  result_df <- as.data.frame(results[[contrast]])

  result_df[] <- lapply(result_df, function(x) {
    if (is.list(x)) {
      vapply(x, function(y) {
        paste(unlist(y), collapse = "; ")
      }, character(1))
    } else {
      x
    }
  })

  write.csv(
    result_df,
    paste0("limma_", contrast, ".csv"),
    row.names = TRUE
  )
}
```

**Problem #3. Choose one pathology group comparison and look up where the probe is in the genome, and the closest transcriptional start site (TSS) gene for the top 5 probes. What does that tell you about the possible pathways that may be deregulated in that pathology group?**

**Problem #4. Take one of the outputs of the linear model and be creative; use ggplot2 to plot the results in different way. Use colours, shapes, or different types of plots.**