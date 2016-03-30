#!/use/bin/env python

import glob
import os
import re
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Extract statistics from STAR log files")
parser.add_argument("-o", "--out",
                    help="Output file",
                    required=True)
parser.add_argument("logs",
                    metavar="LOGFILE",
                    nargs="+",
                    help="Input STAR log files (*.Log.final.out)")
args = parser.parse_args()

log_files = args.logs

# List of stats we want to collect
stat_list = ["Sample", "InputNum", "AverageLen", "MappedNum", "MappedPer",
             "MultiLociNum", "MultiLociPer", "ManyLociNum", "ManyLociPer",
             "UnmappedMismatchPer", "UnmappedShortPer", "UnmappedOtherPer",
             "ChimericNum", "ChimericPer"]

# Create dictionary for storing stats
stats = {}
for stat in stat_list:
    stats[stat] = []

# Map from description in STAR log to stat name
stat_names = {"Number of input reads" : "InputNum",
              "Average input read length" : "AverageLen",
              "Uniquely mapped reads number" : "MappedNum",
              "Uniquely mapped reads %" : "MappedPer",
              "Number of reads mapped to multiple loci" : "MultiLociNum",
              "% of reads mapped to multiple loci" : "MultiLociPer",
              "Number of reads mapped to too many loci" : "ManyLociNum",
              "% of reads mapped to too many loci" : "ManyLociPer",
              "% of reads unmapped: too many mismatches" : "UnmappedMismatchPer",
              "% of reads unmapped: too short" : "UnmappedShortPer",
              "% of reads unmapped: other" : "UnmappedOtherPer",
              "Number of chimeric reads" : "ChimericNum",
              "% of chimeric reads" : "ChimericPer"}

# Extract stats from log files
for log_file in log_files:
    sample = log_file.split("/")[-1].split(".")[0]
    stats["Sample"].append(sample)
    with open(log_file, "r") as log:
        for line in log:
            line = line.strip().split("|")
            if len(line) == 2:
                name = line[0].strip()
                value = line[1].strip().strip("%")
                if name in stat_names.keys():
                    stat = stat_names[name]
                    stats[stat].append(value)

# Convert to Dataframe for easy output
output = pd.DataFrame(stats)
output = output[stat_list]
output.to_csv(args.out, index = False)
