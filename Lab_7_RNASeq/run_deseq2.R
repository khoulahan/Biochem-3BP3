#### RUN DESEQ2 ###################################################################################
# run differential abundance analysis

### PREAMBLE ######################################################################################
library(DESeq2)
library(purrr)
library(dplyr)

date <- Sys.Date()
### READ IN AND REFORMAT COUNTS ###################################################################
read_in_and_reformat_counts <- function() {
    # find feature counts file names
    files <- list.files(
        "feature_counts/", 
        pattern = "feature_counts\\.txt$", 
        full.names = TRUE
        )
    # read in file name
    df_list <- sapply(files, function(x) {
        tmp <- read.delim(x, as.is = TRUE, skip = 1)
        return(tmp[,c(1,7)])
        },
        simplify = FALSE)
    # merge files by gene id
    merged <- reduce(df_list, full_join, by = "Geneid")
    colnames(merged) <- gsub("hisat2\\.|\\.bam", "", colnames(merged))
    return(merged)
    }

### MAIN ##########################################################################################
# read in and reformat gene counts 
gene_counts <- read_in_and_reformat_counts()

# create metadata 
metadata <- data.frame(
    id = colnames(gene_counts)[-1],
    treatment = gsub("HLE_|_1|_2|_3", "", colnames(gene_counts)[-1])
    )

# construct DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = gene_counts,
    colData = metadata,
    design = ~treatment # this tells DESeq how to group the samples for the comparison
    )

# run DESeq2
dds <- DESeq(dds)

# export results
res <- results(dds)
# remove nas
res_nona <- res[!is.na(res$log2FoldChange),]
# order by log2FC
res_nona <- res_nona[order(-abs(res_nona$log2FoldChange)),]

# write to file
write.table(
    res_nona,
    file = paste0(date, '_deseq2_results.txt'),
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
    )

# create normalized counts
vsdata <- vst(dds)
# convert to matrix
vsdata_matrix <- assay(vsdata)

# write to file
write.table(
    vsdata_matrix,
    file = paste0(date, '_normalized_vst_data.txt'),
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
    )

# create pca plot 
png(filename = paste0(date, '_pca_plot.png'), width = 5, height = 5, units = 'in', res = 600)
plotPCA(vsdata, intgroup="treatment")
dev.off()
