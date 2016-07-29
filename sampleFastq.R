#!/usr/bin/env Rscript

# R script for sampling reads from FASTQ files.
# See help for usage.

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("ShortRead"))

parser <- ArgumentParser(description = "Sample reads from FASTQ files")
parser$add_argument("-1", "--pair1",
                    nargs = "+",
                    help = "FASTQ files containing the first read in each pair",
                    required = TRUE)
parser$add_argument("-2", "--pair2",
                    nargs = "*",
                    help = "FASTQ files containing the second read in each pair")
parser$add_argument("-n", "--nreads",
                    type = "integer",
                    default = 1e6,
                    help = "Number of reads to sample")
parser$add_argument("-s", "--seed",
                    type = "integer",
                    default = 1,
                    help = "Random seed for sampling")
parser$add_argument("-o", "--outpath",
                    default = getwd(),
                    help = "Output directory for sampled files")
args <- parser$parse_args()

paired <- !is.null(args$pair2)

if (paired) {
    if (!(length(args$pair2) == length(args$pair1))) {
        stop("Same number of files must be provided for pair1 and pair2")
    }
}

if (!dir.exists(args$outpath)) {
    message(paste("Creating output directory at", args$outpath))
    dir.create(args$outpath)
}

if (paired) {
    pairs <- mapply(c, args$pair1, args$pair2, SIMPLIFY = FALSE)
    message("Found ", length(pairs), " pairs")
    for (pair in pairs) {
        message("Sampling pair: ", basename(pair[1]), " ", basename(pair[2]))
        for (file in pair) {
            fq <- FastqSampler(file, n = args$nreads)

            set.seed(args$seed)
            sampled <- yield(fq)

            out <- strsplit(basename(file), ".", fixed = TRUE)[[1]][1]
            out <- paste0(out, "_sampled.fastq.gz")
            out <- file.path(args$outpath, out)
            message("Writing: ", out)
            writeFastq(sampled, out)
            close(fq)
        }
    }
} else {
    message("Found ", length(args$pair1), " files")
    for (file in args$pair1) {
        message("Sampling file: ", basename(file))
        fq <- FastqSampler(file, n = args$nreads)

        set.seed(args$seed)
        sampled <- yield(fq)

        out <- strsplit(basename(file), ".", fixed = TRUE)[[1]][1]
        out <- paste0(out, "_sampled.fastq.gz")
        out <- file.path(args$outpath, out)
        message("Writing: ", out)
        writeFastq(sampled, out)
        close(fq)
    }
}
message("Done!")
