#!/usr/bin/env python
"""
Run kallisto on multiple files
"""

import argparse
import os
import re
import shlex
import subprocess
import sys

def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Run kallisto on multiple samples",
        epilog="Any additional, unknown, arguments will be passed to kallisto. These must go after the input files.")
    parser.add_argument("inputs",
                        metavar="FASTQ",
                        nargs="+",
                        help="Input FASTQ files")
    parser.add_argument("-i", "--index",
                        help="Path to kallisto index",
                        required=True)
    parser.add_argument("-o", "--output-dir",
                        help="Path to output directory",
                        required=True)
    parser.add_argument("-1", "--regex1",
                        help="Regex for selecting first file in pair",
                        default="_1")
    parser.add_argument("-2", "--regex2",
                        help="Regex for selecting second file in pair",
                        default="_2")
    args, unknown = parser.parse_known_args()
    kal_opts = " ".join(unknown)

    return args, kal_opts


def pair_inputs(inputs, regex1, regex2):
    """
    Take a list of inputs and pair them based on two regular expressions
    """

    inputs1 = []
    inputs2 = []

    for input_file in inputs:
        if re.search(regex1, input_file):
            inputs1.append(input_file)
        elif re.search(regex2, input_file):
            inputs2.append(input_file)
        else:
            sys.exit("""An input did not match either regex. You may need to
                     adjust your patterns""")

    if len(inputs1) == len(inputs2):
        pairs = zip(inputs1, inputs2)
        return pairs
    else:
        sys.exit("""Different numbers of inputs matched each regualar
                 expression. You may need to adjust your patterns.""")


def main():
    """
    Run main code
    """

    args, kal_opts = get_args()

    input_pairs = pair_inputs(args.inputs, args.regex1, args.regex2)

    failed = 0

    for pair in input_pairs:
        sample = os.path.basename(pair[0]))
        sample = sample.split(".")[0]
        print sample
        cmd_str = "kallisto quant --index " + args.index + " --output-dir "
        cmd_str += args.output_dir + "/" + sample + " " + kal_opts + " "
        cmd_str += pair[0] + " " + pair[1]
        cmd = shlex.split(cmd_str)
        failed = subprocess.call(cmd)
        if failed:
            sys.exit("kallisto command '" + cmd_str + "' failed to run")


if __name__ == "__main__":
    main()
