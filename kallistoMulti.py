#!/usr/bin/env python
"""
Run kallisto on multiple files
"""

import argparse
import logging
import os
import shlex
import subprocess
import sys

def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Run kallisto on multiple samples",
        epilog="Any additional, unknown, arguments will be passed to kallisto.")
    #parser.add_argument("inputs",
    #                    metavar="FASTQ",
    #                    nargs="+",
    #                    help="Input FASTQ files")
    parser.add_argument("-i", "--index",
                        help="Path to kallisto index",
                        required=True)
    parser.add_argument("-o", "--output-dir",
                        help="Path to output directory",
                        required=True)
    parser.add_argument("-1", "--read1",
                        help="FASTQ files containing first read in pair",
                        nargs="+")
    parser.add_argument("-2", "--read2",
                        help="FASTQ files containing second read in pair",
                        nargs="+")
    args, unknown = parser.parse_known_args()
    kal_opts = " ".join(unknown)

    return args, kal_opts


def setup_logging(log_file):
    """
    Setup logging system.

    Log is written to 'alignmentStats.log'.
    """

    logger = logging.getLogger("kalMulti")
    logger.setLevel(logging.DEBUG)

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


def main():
    """
    Run main code
    """

    args, kal_opts = get_args()

    log_file = os.path.join(args.output_dir, "kallistoMulti.log")

    setup_logging(log_file)

    logger = logging.getLogger("kalMult." + __name__)

    #input_pairs = pair_inputs(args.inputs, args.regex1, args.regex2)
    input_pairs = zip(args.read1, args.read2)

    failed = 0
    sample_count = 0

    logger.info("Beginning kallisto quantification on multiple samples...\n")
    for pair in input_pairs:
        sample = os.path.basename(pair[0])
        sample = sample.split(".")[0]
        sample_count += 1
        logger.info("Quantifying sample" + sample_count + "of" + len(input_pairs) + "-" + sample)
        cmd_str = "kallisto quant --index " + args.index + " --output-dir "
        cmd_str += args.output_dir + "/" + sample + " " + kal_opts + " "
        cmd_str += pair[0] + " " + pair[1]
        cmd = shlex.split(cmd_str)
        failed = subprocess.call(cmd)
        if failed:
            sys.exit("kallisto command '" + cmd_str + "' failed to run")
        logger.info("Sample" + sample + "complete\n")
    logger.info("Done!")


if __name__ == "__main__":
    main()
