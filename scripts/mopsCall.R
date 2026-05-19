#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    BiocManager::install()
}

if (!require("cn.mops", quietly = TRUE)){
    BiocManager::install("cn.mops")
}

if(!require("magrittr", quietly=TRUE)){
    install.packages("magrittr", repos = "http://cran.us.r-project.org")
}

suppressMessages(library(cn.mops))
suppressMessages(library(magrittr))

options(scipen = 999)

# Get arguments.
args <- commandArgs(trailingOnly = TRUE)
result_dir <- args[1]
bin_size <- as.integer(args[2])
threads <- as.integer(args[3])
bam_files <- tail(args, -3)

# Get read counts from BAM files (all chromosomes, no species-specific filtering)
bam_data_ranges <- getReadCountsFromBAM(bam_files, WL = bin_size, parallel = threads)

# Remove chromosomes where any sample has zero average read count
# (cn.MOPS normalizeChromosomes requires non-zero average per chr per sample)
count_matrix <- as.matrix(mcols(bam_data_ranges))
chr_names <- as.character(seqnames(bam_data_ranges))
chr_sample_means <- sapply(split(seq_along(chr_names), chr_names), function(idx) {
    colMeans(count_matrix[idx, , drop = FALSE])
})
# chr_sample_means: rows=samples, cols=chromosomes
good_chrs <- colnames(chr_sample_means)[apply(chr_sample_means, 2, function(x) all(x > 0))]
cat("Keeping chromosomes:", paste(good_chrs, collapse = ", "), "\n")
bam_data_ranges <- bam_data_ranges[as.character(seqnames(bam_data_ranges)) %in% good_chrs, ]

# Call CNVs using haplocn.mops for haploid genomes
results <- haplocn.mops(bam_data_ranges, parallel = threads) %>% calcIntegerCopyNumbers()
cnvs <- cnvs(results)

# Format GRanges to bed-like data.frame
granges_to_bed <- function(gr) {
    bed <- data.frame(
        chrom = as.character(seqnames(gr)),
        chromStart = start(ranges(gr)),
        chromEnd = end(ranges(gr)),
        copyNumber = sub("^CN", "", gr$CN),
        median = gr$median
    )
    return(bed)
}

write_bed <- function(bed, file_name) {
    write.table(bed, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_sample_cnvs <- function(sample_name, cnv_results) {
    sample_cnvs <- cnv_results[cnv_results$sampleName == sample_name]
    bed <- granges_to_bed(sample_cnvs)
    file_name <- sub(".bam$", ".temp.bed", sample_name)
    path <- paste0(result_dir, file_name)
    write_bed(bed, path)
}

# Output CNV bed files for all input samples.
sapply(basename(bam_files), write_sample_cnvs, cnv_results = cnvs)
