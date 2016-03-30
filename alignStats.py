#!/usr/bin/env python
"""
Extract alignment statistics from a SAM/BAM file.

Adapted from the Celloline stats script
available at: https://github.com/Teichlab/celloline/blob/master/lib/stats.py
"""

import os
import sys
import re
import argparse
import pysam
import logging
import cPickle as pickle
from collections import Counter, defaultdict, OrderedDict
from intervaltree import IntervalTree
from joblib import Parallel, delayed

#LOAD GTF FILE
def load_gtf(gtf_path):
    """
    Load a GTF annotation and create an index using IntervalTrees.

    Args:
        gtf_path: Path to the GTF file to load.

    Returns:
        Dictionary containing IntervalTree indexes of the annotation.
    """

    gtf_index = defaultdict()
    with open(gtf_path) as gtf_file:
        for line in gtf_file:
            if not line.startswith("#"):
                entry = line.split("\t")
                entry_addition = entry[8]
                entry_addition = entry_addition.split(";")
                entry_addition = entry_addition[0].split(" ")
                gene_id = entry_addition[1]

                feature = entry[2]
                #TYPE(Gene, exon etc.), START, END, STRAND, gene_ID
                info = [feature, entry[3], entry[4], entry[6], gene_id]

                #Build GTF INDEX
                if feature != "" and entry[3] != entry[4]:
                    if entry[0] in gtf_index:
                        index = gtf_index[entry[0]]
                    else:
                        index = IntervalTree()
                    index.addi(int(info[1]), int(info[2]), info)
                    gtf_index[entry[0]] = index

    return gtf_index


def gen_stats(input_file, input_type, sample_name, gtf_dict):
    """
    Generate alignment stats from a SAM/BAM file.

    Loop over alignments in a SAM/BAM file and extract statistics such as the
    numer of reads aligned to introns, exons, intergenic regions etc.

    Args:
        input_file: An open BAM or SAM file.
        input_type: Whether the file is 'bam' or 'sam'.
        sample_name: A name relating to this file.
        gtf_dict: Dictionary containing GTF index.

    Returns:
        Dictionary containing alignment statistics.
    """

    logger = logging.getLogger("stats." + sample_name[0:10])

    #OUTPUT TABLE CONTAING STATS
    output_table = OrderedDict()

    #Dict indicating to which genes a specific read maps to
    #It is a temporary dict
    exonic_mappings_temp = defaultdict(str)

    #Dict indicating which read is multi-mapped
    #It is a temporary dict
    exonic_multi_table = defaultdict(str)

    # Sample
    output_table["sample"] = sample_name

    #MAPPABILITY
    output_table["total"] = 0
    output_table["mapped"] = 0
    output_table["unmapped"] = 0
    output_table["unique"] = 0
    output_table["multi"] = 0

    #CODING VERSUS NON-CODING REGIONS
    output_table["intergenic"] = 0
    output_table["intragenic"] = 0
    output_table["exonic"] = 0
    output_table["intronic"] = 0
    output_table["ambigious"] = 0

    #CODING REGIONS MAPPABILITY
    output_table["exonicU"] = 0
    output_table["exonicM"] = 0

    #ALIGNMENT CODING VS NONCODING
    output_table["alignments"] = 0
    output_table["multi-intergenic"] = 0
    output_table["multi-intragenic"] = 0
    output_table["multi-exonic"] = 0
    output_table["multi-intronic"] = 0
    output_table["multi-ambigious"] = 0

    #ERROR
    output_table["perfect"] = 0
    output_table["partly_perfect"] = 0
    output_table["mapped_no_correct"] = 0
    for i in range(0, 10):
        output_table["S_" + str(i)] = 0
    output_table["S_10+"] = 0
    output_table["I"] = 0
    output_table["D"] = 0
    output_table["INDEL"] = 0

    reads = Counter()

    if input_type == "bam":
        ref_map = input_file.references
        input_file = input_file.fetch(until_eof=True)

    line_count = 0

    for line in input_file:

        line_count += 1

        if input_type == "bam":                 # BAM input line
            split = str(line).split("\t")
            split[2] = ref_map[int(split[2])]
            split[3] = int(split[3]) + 1
        elif not line.startswith("@"):         # SAM input line
            split = line.split("\t")
        else:
            continue

        read_name = split[0]
        flag_code = int(split[1])
        chrom = split[2]
        pos = split[3]
        errors = split[5]

        errors_a = list(errors)
        number = ""
        num = 0
        error_table = defaultdict(int)
        name_and_flag = read_name

        #CHECK IF READ MAPPED OR UNMAPPED
        #IT US UNMAPPED
        if flag_code & 0x0004 != 0:
            output_table["unmapped"] += 1
            output_table["total"] += 1
            error_table["*"] += 1
        #IT IS MAPPED
        else:
            if flag_code & 0x0001 != 0:         #This is paired end
                if flag_code & 0x0040 != 0:     #1st read
                    name_and_flag += ";first"
                if flag_code & 0x0080 != 0:     #2nd read
                    name_and_flag += ";second"

            # CHECK TO WHICH GENE(S) IT MAPPED TO
            genes_info, num_genes, num_exons = get_gene(gtf_dict, [chrom, pos])

            output_table["alignments"] += 1.0

            #STATS
            if name_and_flag not in reads:
                reads[name_and_flag] += 1
                output_table["unique"] += 1
                output_table["total"] += 1
                output_table["mapped"] += 1

                if num_genes == 0:
                    output_table["intergenic"] += 1
                elif num_genes == 1:
                    output_table["intragenic"] += 1
                    if num_exons == 0:
                        output_table["intronic"] += 1
                    else:
                        output_table["exonic"] += 1
                        output_table["exonicU"] += 1
                        exons = []
                        if name_and_flag in exonic_mappings_temp:
                            exons = exonic_mappings_temp[name_and_flag]
                        exons.append([genes_info[0], chrom, pos])
                        exonic_mappings_temp[name_and_flag] = exons
                elif num_genes > 1:
                    output_table["ambigious"] += 1

            #READ IS MULTI-MAPPED
            else:
                if reads[name_and_flag] == 1:
                    output_table["unique"] -= 1
                    output_table["exonicU"] -= 1
                    output_table["multi"] += 1
                reads[name_and_flag] += 1
                exons = []

                #GET KNOWLEDGE IF FIRST MAPPING EXONIC OR INTRONIC
                if name_and_flag in exonic_mappings_temp:
                    exons = exonic_mappings_temp[name_and_flag]
                if num_genes == 0:
                    output_table["multi-intergenic"] += (1)
                elif num_genes == 1:
                    output_table["multi-intragenic"] += (1)
                    if num_exons == 0:
                        output_table["multi-intronic"] += (1)
                    else:
                        output_table["multi-exonic"] += (1)
                        exons.append([genes_info[0], chrom, pos])
                elif num_genes > 1:
                    output_table["multi-ambigious"] += (1)
                #IF AT LEAST ONE EXONIC ALIGNMENT
                if len(exons) > 0:
                    exonic_multi_table[name_and_flag] = exons

            #PARSE MAPPING ERRORS
            for i in errors_a:
                if re.match("[0-9]", i):
                    number += (i)
                elif re.match("[A-Z]", i):
                    num = int(number)
                    error_table[i] += num
                    number = ""

            #TABLE OF HOW MANY READS MAP PERFECT, PARTLY PERFECT ETC
            if "M" in  error_table and len(error_table) == 1:
                output_table["perfect"] += 1
            elif "M" in error_table and len(error_table) > 1:
                output_table["partly_perfect"] += 1
            elif "M" not in error_table and "*" not in error_table:
                output_table["mapped_no_correct"] += 1

            if "S" in error_table:
                if int(error_table["S"]) < 10:
                    output_table["S_" + str(error_table["S"])] += 1
                else:
                    output_table["S_10+"] += 1
            elif "S" not in error_table:
                output_table["S_0"] += 1

            if "I" in error_table:
                output_table["I"] += 1

            if "D" in error_table:
                output_table["D"] += 1

            if "I" in error_table or "D" in error_table:
                output_table["INDEL"] += 1

        if (line_count % 1000000) == 0:
            logger.debug(sample_name + " line " + str(line_count) + "...")

    output_table["exonicM"] = len(exonic_multi_table.keys())

    return output_table


def get_stats_line(stats_table):
    """
    Get an output line from a stats table.

    Args:
        stats_table: Dictionary of alignment statistics.

    Returns:
        String representing the results for one file.
    """

    logger = logging.getLogger("stats.extract")

    out_line = ""

    for stat, value in stats_table.iteritems():
        if stat in ["unique", "multi", "intragenic", "intergenic",
                    "exonic", "intronic", "ambigious", "exonicM", "exonicU"]:
            value = (value + 0.0) / (stats_table["mapped"] + 0.0)
            value = "%.2f" % (100.0 * (value))
        elif stat in ["multi-intragenic", "multi-intergenic", "multi-exonic",
                      "multi-intronic", "multi-ambigious"]:
            value = (value + 0.0)
            if stats_table["alignments"] != 0:
                value = value / (stats_table["alignments"] + 0.0)
            value = "%.2f" % (100.0 * (value))

        value = str(value)

        if not stat == "sample":
            out_line += "," + value
        else:
            out_line += value

        logger.debug(stat + " : " + value)

    out_line += "\n"

    return out_line


def write_stats(output_path, stats_list):
    """
    Write a series of results to a file.

    Args:
        output_path: Path to write results to.
        stats_list: List of dictionaries containing results from input files.
    """

    cols = stats_list[0].keys()

    with open(output_path, "w") as out_file:
        out_file.write(",".join(cols) + "\n")
        for stats_table in stats_list:
            stats_line = get_stats_line(stats_table)
            out_file.write(stats_line)


def get_gene(gtf_dict, pos_pair):
    """
    Identify which genes overlap a given position.

    Args:
        gtf_dict: Dictionary containing GTF index.
        pos_pair: Tuple containing genomic position (chrom, pos).

    Returns:
        Tuple containing the list of overlapping genes, the number of
        overlapping genes and the number of overlapping exons.
    """

    num_genes = 0
    num_exons = 0

    if pos_pair[0] not in gtf_dict:
        #print ("Ignored pos: " + pos_pair[0])
        return ([], num_genes, num_exons)

    entries = gtf_dict[pos_pair[0]]
    pos = int(pos_pair[1])

    found = []
    found = entries.search(pos)

    gene_list = []
    for entry in found:
        info = entry[2]
        if info[0] == "gene":
            gene_list.append(info)
            num_genes += 1
        elif info[0] == "exon":
            num_exons += 1

    return (gene_list, num_genes, num_exons)


def process_file(input_file, input_type, index, is_parallel):
    """
    Process an individual SAM/BAM file.

    How we want to process the file depends on the input type and whether we
    are operating in parallel. If in parallel the index must be loaded for each
    input file. If the input is a BAM file it needs to be read using Pysam, if
    SAM it can be read directly as a text file.

    Args:
        input_file: Path to the input file.
        input_type: Whether the file is 'bam' or 'sam'.
        index: If operating in parallel a string to the index file, if not the
               loaded GTF index dictionary.
        is_parallel: Whether to operate in parallel.

    Returns:
        Dictionary containing alignment statistics for the input file.
    """

    sample_name = input_file.split("/")[-1]

    logger = logging.getLogger("stats." + sample_name[0:10])

    logger.info("Processing " + sample_name + "...")

    if is_parallel:
        logger.info("Loading index...")
        with open(index, "rb") as index_file:
            loaded_index = pickle.load(index_file)
        logger.info("Loaded.")
    else:
        loaded_index = index

    if input_type == "sam":
        logger.info("Parsing SAM file...")
        with open(input_file) as sam:
            output_table = gen_stats(sam, input_type, sample_name, loaded_index)
    elif input_type == "bam":
        logger.info("Parsing BAM file...")
        bam = pysam.AlignmentFile(input_file, "rb")
        output_table = gen_stats(bam, input_type, sample_name, loaded_index)

    logger.info("Finished " + sample_name)

    return output_table


def get_index(args):
    """
    Load a GTF index if available or create from GTF file if not found.

    If a valid path to an index file is given that file will be loaded. If no
    index file was specified or the file does not exist the annotation will be
    read from a GTF file. It will then be pickled if an index file is specified.

    When running in parallel the path to the index file is returned rather than
    the index dictionary itself.

    Args:
        args: Options from the command line.

    Returns:
        Dictionary containing GTF index or path to index file if in parallel.
    """

    logger = logging.getLogger("stats.index")

    if args.index and os.path.isfile(args.index):
        logger.info("Index found at " + args.index)

        if not args.is_parallel:
            logger.info("Loading index...")
            with open(args.index, "rb") as index_file:
                index = pickle.load(index_file)
            logger.info("Loaded.")
        else:
            index = args.index

    elif args.gtf and os.path.isfile(args.gtf):

        logger.info("No index file found.")
        logger.info("Loading GTF file...")
        gtf_dict = load_gtf(args.gtf)
        logger.info("Loaded.")

        if args.index:
            logger.info("Saving index to " + args.index + "...")
            with open(args.index, "wb") as index_file:
                pickle.dump(gtf_dict, index_file, -1)
            logger.info("Saved.")

        if not args.is_parallel:
            index = gtf_dict
        else:
            index = args.index

    return index


def get_args():
    """
    Read arguments from the command line and check they are valid.
    """

    logger = logging.getLogger("stats.args")

    parser = argparse.ArgumentParser(
        description="Extract alignment statistics from a SAM/BAM file")
    parser.add_argument("inputs",
                        metavar="SAM/BAM",
                        nargs="+",
                        help="Input SAM or BAM files")
    parser.add_argument("-o", "--out",
                        help="Output file",
                        required=True)
    parser.add_argument("-g", "--gtf",
                        help="GTF annotation file")
    parser.add_argument("-i", "--index",
                        help="""Annotation index file. Required when
                                operating in parallel.""")
    parser.add_argument("-t", "--type",
                        choices=["sam", "bam"],
                        help="Type of input file",
                        required=True)
    parser.add_argument("-p", "--parallel",
                        type=int,
                        default=1,
                        help="""Number of files to process in parallel.
                             Requires N + 1 threads if greater than 1.""")
    args = parser.parse_args()

    args.is_parallel = False

    if args.parallel < 1:
        logger.error("Number of parallel files must be positive")
        sys.exit()
    elif args.parallel > 1:
        args.is_parallel = True
        logger.info("Running with " + str(args.parallel) + " jobs")

    if args.is_parallel and not args.index:
        logger.error("Index file is required when running in parallel.")
        sys.exit()

    if not (args.index and os.path.isfile(args.index)):
        if not (args.gtf and os.path.isfile(args.gtf)):
            logger.error("No GTF or index file found.")
            sys.exit()

    return args


def setup_logging():
    """
    Setup logging system.

    Log is written to 'alignmentStats.log'.
    """

    logger = logging.getLogger("stats")
    logger.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    file_handler = logging.FileHandler('alignmentStats.log')
    file_handler.setLevel(logging.INFO)

    # create console handler with a higher log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    format_str = "[%(asctime)s] %(levelname)s %(name)s: %(message)s"
    formatter = logging.Formatter(format_str, "%Y-%m-%d %H:%M:%S")
    file_handler.setFormatter(formatter)
    format_str = "[%(asctime)s] %(message)s"
    formatter = logging.Formatter(format_str, "%H:%M:%S")
    console_handler.setFormatter(formatter)

    # add the handlers to logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)


def main():
    """
    Main function.

    1. Setup logging
    2. Get arguments
    3. Get index
    4. Process files
    5. Write output
    """

    setup_logging()

    logger = logging.getLogger("stats." + __name__)

    args = get_args()

    index = get_index(args)

    logger.warning("Positions not in annotation will be ignored.")

    logger.info("Found " + str(len(args.inputs)) + " input file(s):")
    for input_file in sorted(args.inputs):
        logger.debug(input_file)

    if args.is_parallel:
        stats = Parallel(n_jobs=args.parallel,
                         verbose=100)(delayed(process_file)(input_file,
                                                            args.type,
                                                            index,
                                                            args.is_parallel)
                                      for input_file in args.inputs)
    else:
        stats = []
        for input_file in args.inputs:
            output_table = process_file(input_file, args.type, index,
                                        args.is_parallel)
            stats.append(output_table)

    write_stats(args.out, stats)

if __name__ == "__main__":
    main()
