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


def setup_logging(output_dir):
    """
    Setup logging system.

    Log is written to 'alignmentStats.log'.
    """

    logger = logging.getLogger("kalMult")
    logger.setLevel(logging.DEBUG)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    log_file = os.path.join(output_dir, "kallistoMulti.log")

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

    setup_logging(args.output_dir)

    logger = logging.getLogger("kalMult." + __name__)

    input_pairs = zip(args.read1, args.read2)

    failed = 0
    sample_count = 0

    logger.info("Beginning kallisto quantification on multiple samples...")
    for pair in input_pairs:
        sample = os.path.basename(pair[0])
        sample = sample.split(".")[0]
        sample_count += 1
        logger.info("Quantifying sample " + sample + " (" + str(sample_count) + " of " +
                str(len(input_pairs)) + ")...")

        cmd_str = "kallisto quant --index " + args.index + " --output-dir "
        cmd_str += args.output_dir + "/" + sample + " " + kal_opts + " "
        cmd_str += pair[0] + " " + pair[1]
        cmd = shlex.split(cmd_str)

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        kal_out = ""
        while True:
            out = process.stderr.read(1)
            if out == '' and process.poll() != None:
                break
            if out != '':
                kal_out += out
                sys.stdout.write(out)
                sys.stdout.flush()

        kal_log_path = os.path.join(args.output_dir, sample, "kallisto.log")
        with open(kal_log_path, "w") as kal_log:
            logger.info("Writing kallisto log to " + kal_log_path + "...")
            kal_log.write(kal_out)

        logger.info("Sample " + sample + " complete")
    logger.info("Done!")


if __name__ == "__main__":
    main()
