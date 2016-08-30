#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description="Merge FASTQ files")
parser.add_argument("-s", "--separator",
                    help="Filename separator",
                    default="_")
parser.add_argument("-m", "--merge",
                    help="Section of the filename to merge",
                    type=int,
                    default=1)
parser.add_argument("-o", "--outdir",
                    help="Path to output directory for merged files",
                    default=os.getcwd())
parser.add_argument("fastqs",
                    metavar="FASTQ",
                    nargs="+",
                    help="FASTQ files to merge")
args = parser.parse_args()

print(args.fastqs)
