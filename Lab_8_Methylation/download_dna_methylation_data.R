#### DOWNLOAD DNA METHYLATION DATA ################################################################
# download dna methylation data from GEO
# accession: GSE100197

### PREAMBLE ######################################################################################
library(rlang) 
library(GEOquery)

date <- Sys.Date()
### DOWNLOAD METHYLATION DATA ###################################################################
download_methylation_data <- function() {
    # Download and read the GEO series matrix for GSE100197
    url <- paste0(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/",
        "GSE100nnn/GSE100197/matrix/",
        "GSE100197_series_matrix.txt.gz"
        )

    filename <- "GSE100197_series_matrix.txt.gz"

    if (!file.exists(filename)) {
    download.file(
        url,
        destfile = filename,
        mode = "wb",
        method = "libcurl"
    )
    }

    #Check that the file exists and contains data
    stopifnot(file.exists(filename), file.info(filename)$size > 0)

    # Read the downloaded compressed file
    con <- gzfile(filename, open = "rt")
    lines <- readLines(con, warn = FALSE)
    close(con)

    # Extract the data table
    start <- grep("^!series_matrix_table_begin", lines)
    end <- grep("^!series_matrix_table_end", lines)

    if (length(start) != 1 || length(end) != 1 || end <= start + 1) {
    stop("No valid data table found. Check the downloaded GEO series matrix.")
    }

    # Write only the table to file
    output_filename <- 'GSE100197_methylation_data.txt'
    writeLines(lines[(start + 1):(end - 1)], output_filename)
    }

### DOWNLOAD METHYLATION METADATA #################################################################
download_methylation_metadata <- function() {
    # Download and extract sample-level clinical metadata from GEO
    gse <- "GSE100197"
    geo_dir <- getwd()
    soft_file <- paste0(gse, "_family.soft.gz")

    # Download directly from NCBI GEO
    geo_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    "GSE100nnn/", gse, "/soft/",
    gse, "_family.soft.gz"
    )

    if (!file.exists(soft_file)) {
    download.file(
        geo_url,
        destfile = soft_file,
        mode = "wb",
        method = "libcurl"
        )
    }

    stopifnot(file.exists(soft_file), file.info(soft_file)$size > 0)

    # Read the downloaded GEO sample records
    con <- gzfile(soft_file, open = "rt")
    soft_lines <- readLines(con, warn = FALSE)
    close(con)

    # return downloaded file 
    return(soft_lines)
    }

### REFORMAT METHYLATION METADATA #################################################################
reformat_methylation_metadata <- function(soft_lines) {
    gse <- "GSE100197"
    # Split the SOFT file into individual GSM sample records
    sample_starts <- grep("^\\^SAMPLE = GSM", soft_lines)

    if (length(sample_starts) == 0) {
    stop("No GSM sample records found in the downloaded GEO file.")
    }

    sample_ends <- c(sample_starts[-1] - 1L, length(soft_lines))

    # Extract sample-level metadata and clinical characteristics
    metadata_list <- lapply(seq_along(sample_starts), function(i) {

    block <- soft_lines[sample_starts[i]:sample_ends[i]]

    # Extract a GEO field, retaining repeated values
    get_field <- function(field) {
        pattern <- paste0("^!Sample_", field, " = ")
        x <- grep(pattern, block, value = TRUE)
        sub(pattern, "", x)
    }

    accession <- get_field("geo_accession")
    title <- get_field("title")
    source <- get_field("source_name_ch1")
    characteristics <- get_field("characteristics_ch1")

    # Store basic sample information
    result <- list(
        geo_accession = if (length(accession)) accession[1] else NA_character_,
        title = if (length(title)) title[1] else NA_character_,
        source_name = if (length(source)) source[1] else NA_character_,
        characteristics_raw = paste(characteristics, collapse = " | ")
    )

    # Convert GEO characteristic entries into separate clinical columns
    for (entry in characteristics) {

        parts <- strsplit(entry, ":", fixed = TRUE)[[1]]

        if (length(parts) >= 2) {
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = ":"))

        # Make column names consistent
        key <- make.names(key)

        if (key %in% names(result)) {
            result[[key]] <- paste(
            result[[key]], value, sep = "; "
            )
        } else {
            result[[key]] <- value
        }
        }
    }

    result
    })

    # Combine records, preserving all available clinical fields
    all_fields <- unique(unlist(lapply(metadata_list, names)))

    # Create the metadata data.frame 
    metadata <- as.data.frame(
    do.call(rbind, lapply(metadata_list, function(x) {
        row <- setNames(
        as.list(rep(NA_character_, length(all_fields))),
        all_fields
        )

        row[names(x)] <- lapply(x, function(y) {
        paste(as.character(unlist(y)), collapse = "; ")
        })

        row
    })),
    stringsAsFactors = FALSE,
    check.names = FALSE
    )

    # Convert any remaining list columns to character strings
    metadata[] <- lapply(metadata, function(x) {
    if (is.list(x)) {
        vapply(x, function(y) {
        paste(as.character(unlist(y)), collapse = "; ")
        }, character(1))
    } else {
        as.character(x)
    }
    })

    # Save the clinical metadata
    write.csv(
        metadata,
        paste0(gse, "_clinical_metadata.csv"),
        row.names = FALSE,
        na = ""
        )
}
### MAIN ##########################################################################################
# download methylation data
download_methylation_data()

# download methylation metadata
soft_lines <- download_methylation_metadata()

# reformat methylation metadata
reformat_methylation_metadata(soft_lines)