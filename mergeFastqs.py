#!/usr/bin/env python
"""
Merge multiple FASTQ files based on a filename pattern.
"""

import argparse
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
                             "K = Keep this section, " +
                             "M = Merge using this section, " +
                             "R = Remove this section. " +
                             "For example if the filname structure was: " +
                             "'SAMPLE_READ_LANE_DATE', " +
                             "to merge on LANE and remove DATE the pattern " +
                             "would be 'K_K_M_R'.",
                        required=True)
    parser.add_argument("fastqs",
                        metavar="FASTQ",
                        nargs="+",
                        help="FASTQ files to merge")
    args = parser.parse_args()

    return args


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


def process_pattern(pat, sep):
    """
    Process the pattern argument to determine what to do with each filename
    section.
    """

    print("Input pattern is: " + pat)
    split_pat = pat.split(sep)
    print("Split pattern is: " + str(split_pat))
    merge_pat = []
    for section in split_pat:
        if (section == "K" or section == "M"):
            merge_pat.append(section)
    print("Merged pattern is: " + str(merge_pat))

    return split_pat


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

    for groupname, filenames in groups.iteritems():
        print("Merging group " + groupname + "...")
        outpath = os.path.join(outdir, groupname)
        print("Creating merge file " + outpath + "...")
        with open(outpath, "wb") as outfile:
            for filename in filenames:
                print("Adding file " + filename + "...")
                with open(filename, "rb") as fq_file:
                    shutil.copyfileobj(fq_file, outfile)


def main():
    """
    Run main code
    """

    args = get_args()
    pattern = process_pattern(args.pattern, args.separator)

    file_groups = group_filenames(args.fastqs, pattern, args.separator)

    merge_files(file_groups, args.outdir)


if __name__ == "__main__":
    main()
