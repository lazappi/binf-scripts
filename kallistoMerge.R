#!/usr/bin/env Rscript

# R script for merging kallisto output.
# See help for usage.

suppressPackageStartupMessages(library("argparse"))
library(tximport)
library(readr)

parser <- ArgumentParser(description = "Merge kallisto output")
parser$add_argument("-d", "--directory",
                    required = TRUE,
                    help = "Directory to search for kallisto output")
parser$add_argument("-t", "--tx2gene",
                    required = TRUE,
                    help = "Mapping from transcripts to genes (TSV)")
parser$add_argument("-o", "--outpath",
                    default = getwd(),
                    help = "Output directory")
args <- parser$parse_args()

message("Searching ", args$directory, " for abundance.tsv files...")
tsv.files <- list.files(args$directory, "abundance.tsv", recursive = TRUE)

if (length(tsv.files) > 0) {
    message("Found ", length(tsv.files), " abundance.tsv files")
} else {
    stop("No abundance.tsv files found")
}

message("Reading tx2gene mapping...")
tx2gene <- read.delim(args$tx2gene)

message("Reading transcript abundances...")
txi.tx <- tximport(file.path(args$directory, tsv.files), type = "kallisto",
                   tx2gene = tx2gene, reader = read_tsv, txOut = TRUE)
message("Summarising to genes...")
txi.sum <- summarizeToGene(txi.tx, tx2gene)
message("Summarising to genes (scaledTPM)...")
txi.sum.scl <- summarizeToGene(txi.tx, tx2gene,
                               countsFromAbundance = "scaledTPM")
message("Summarising to genes (lengthScaledTPM)...")
txi.sum.len <- summarizeToGene(txi.tx, tx2gene,
                               countsFromAbundance = "lengthScaledTPM")

if (!dir.exists(args$outpath)) {
    message(paste("Creating output directory at", args$outpath))
    dir.create(args$outpath)
}

message("Writing transcript tables...")

for (name in names(txi.tx)) {
    item <- txi.tx[[name]]
    if (is.matrix(item)) {
        colnames(item) <- dirname(tsv.files)
        out <- paste("kallisto", "tx", name, sep = "_")
        out <- paste0(out, ".tsv")
        out <- file.path(args$outpath, out)
        write.table(item, out, sep = "\t", quote = FALSE)
    }
}

message("Writing gene tables...")

for (name in names(txi.sum)) {
    item <- txi.sum[[name]]
    if (is.matrix(item)) {
        colnames(item) <- dirname(tsv.files)
        out <- paste("kallisto", "gene", name, sep = "_")
        out <- paste0(out, ".tsv")
        out <- file.path(args$outpath, out)
        write.table(item, out, sep = "\t", quote = FALSE)
    }
}

item <- txi.sum.scl$counts
colnames(item) <- dirname(tsv.files)
out <- paste("kallisto", "gene", "scaledTPM", "counts", sep = "_")
out <- paste0(out, ".tsv")
out <- file.path(args$outpath, out)
write.table(item, out, sep = "\t", quote = FALSE)

item <- txi.sum.len$counts
colnames(item) <- dirname(tsv.files)
out <- paste("kallisto", "gene", "lengthScaledTPM", "counts", sep = "_")
out <- paste0(out, ".tsv")
out <- file.path(args$outpath, out)
write.table(item, out, sep = "\t", quote = FALSE)

message("Done!")
