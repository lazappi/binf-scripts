#!/usr/bin/env python
"""
Merge multiple FASTQ files based on a filename pattern.
"""

import argparse
import os

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
                             "Should have the format [CODE][SEP]...[CODE] " +
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


def main():
    """
    Run main code
    """

    args = get_args()


if __name__ == "__main__":
    main()
