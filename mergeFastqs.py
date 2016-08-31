#!/usr/bin/env python
"""
Merge multiple FASTQ files based on a filename pattern.
"""

import argparse
import logging
import os
import sys
import shutil
from collections import defaultdict

def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(description="Merge FASTQ files based on name")
    parser.add_argument("-s", "--separator",
                        help="Filename separator. Default is '_'.",
                        default="_")
    parser.add_argument("-o", "--outdir",
                        help="Path to output directory for merged files.",
                        default=os.getcwd())
    parser.add_argument("-p", "--pattern",
                        help="Pattern used to decide which files to merge." +
                             "Should have the format [CODE][SEP]...[CODE], " +
                             "where SEP is the separator and CODE is one of: " +
                             "K = Keep this section or " +
                             "M = Merge this section, "
                             "For example if the filename structure was: " +
                             "'SAMPLE_READ_LANE_DATE.fastq', " +
                             "to merge on LANE and DATE the pattern " +
                             "would be 'K_K_M_M', which would produce a" +
                             "merged file named 'SAMPLE_DATE.fastq'",
                        required=True)
    parser.add_argument("fastqs",
                        metavar="FASTQ",
                        nargs="+",
                        help="FASTQ files to merge")
    args = parser.parse_args()

    return args


def setup_logging(outdir):
    """
    Setup logging system.

    Log is written to 'mergeFastqs.log'.
    """

    logger = logging.getLogger("mergeFQs")
    logger.setLevel(logging.DEBUG)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    log_file = os.path.join(outdir, "mergeFastqs.log")

    # create file handler which logs even debug messages
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)

    # create console handler with a higher log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

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


def merge_filename(filename, pat, sep):
    """
    Apply a merging pattern to a filename.
    """

    filebase = os.path.basename(filename).split(".")[0]
    ext = ".".join(os.path.basename(filename).split(".")[1:])

    split_file = filebase.split(sep)

    if not len(split_file) == len(pat):
        sys.exit("File " + filename + " does not match pattern " + str(pat))

    merge_file = []
    for idx in xrange(0, len(split_file)):
        code = pat[idx]
        file_sec = split_file[idx]
        if code == "K":
            merge_file.append(file_sec)

    return sep.join(merge_file) + "." + ext


def group_filenames(filenames, pat, sep):
    """
    Group files based on their merged file names.
    """

    groups = defaultdict(list)

    for filename in filenames:
        group = merge_filename(filename, pat, sep)
        groups[group].append(filename)

    return(groups)


def merge_files(groups, outdir):
    """
    Merge files that belong to the same filename group.
    """

    logger = logging.getLogger("mergeFQs." + "merge")

    for groupname, filenames in groups.iteritems():
        logger.info("Merging group " + groupname + " with " + len(filenames) +
                    "files...")
        outpath = os.path.join(outdir, groupname)
        logger.info("Creating merge file " + outpath + "...")
        with open(outpath, "wb") as outfile:
            for filename in filenames:
                logger.info("Adding file " + filename + "...")
                with open(filename, "rb") as fq_file:
                    shutil.copyfileobj(fq_file, outfile)


def main():
    """
    Run main code
    """

    args = get_args()

    setup_logging(args.outdir)
    logger = logging.getLogger("mergeFQs." + __name__)

    logger.info(str(len(args.fastqs)) + " input files provided")
    logger.info("Filename pattern is " + args.pattern)
    pattern = args.pattern.split(args.sep)
    ex_file = args.fastqs[0]
    ex_merge = merge_filename(ex_file, pattern, args.separator)
    logger.info("Example merge: " + ex_file + " -> " +
                os.path.join(args.outdir, ex_merge))
    file_groups = group_filenames(args.fastqs, pattern, args.separator)
    logger.info(len(file_groups) + " file groups found...")
    merge_files(file_groups, args.outdir)


if __name__ == "__main__":
    main()
