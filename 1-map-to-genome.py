#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# A simple
#
#
# Copyright (c) 2014-2016 Nick Semenkovich <semenko@alum.mit.edu>.
#   https://nick.semenkovich.com/
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
#   http://gordonlab.wustl.edu/
#
# This software is released under the MIT License:
#  http://opensource.org/licenses/MIT
#
# Source: https://github.com/GordonLab/riesling-pipeline

from __future__ import absolute_import, division, print_function, unicode_literals

__author__ = 'Nick Semenkovich <semenko@alum.mit.edu>'
__copyright__ = 'Gordon Lab at Washington University in St. Louis'
__license__ = 'MIT'
__version__ = '1.0.3'

from collections import OrderedDict
import _logshim
import _script_helpers
import argparse
import glob
import os
import pprint
import re


def find_paired_ends(input_path, verbose=False):
    """
    Given an input path, return

    :param input_path:
    :return:
    """
    find_pe_logger = _logshim.getLogger('find_paired_ends')

    # TODO: Modularize all this!

    if not os.path.isdir(input_path):
        raise ValueError("Input must be a directory. You gave: %s" % (input_path))

    all_files = glob.glob(input_path + "/*.PE1.fastq.gz")  # Must have .PEX. in title
    all_files.extend(glob.glob(input_path + "/*.PE2.fastq.gz"))
    all_files.extend(glob.glob(input_path + "/*.PE1.fastq"))
    all_files.extend(glob.glob(input_path + "/*.PE2.fastq"))

    if len(all_files) == 0:
        raise ValueError("Input directory is empty!")


    # Given paired ends, we must always have an even number of input files.
    if len(all_files) % 2 != 0:
        raise ValueError("Input directory contains an odd number of files.")

    re_pattern = re.compile(r'^(.*)\.PE(\d)(\.fastq|\.fastq\.gz)$')

    file_dict = OrderedDict()

    prefixes_seen = []
    pe_seen = []
    for file in sorted(all_files):
        if not os.access(file, os.R_OK):
            raise OSError("Cannot read file: %s" % (file))

        filename_only = file.rsplit('/', 1)[-1]
        result = re.match(re_pattern, filename_only)

        file_dict[file] = {'prefix': str(result.group(1)),
                           'PE': int(result.group(2))}

        prefixes_seen.append(file_dict[file]['prefix'])
        pe_seen.append(file_dict[file]['PE'])

    if len(set(pe_seen)) != 2:
        raise ValueError("Saw %d paired ends, expecting exactly two. That's confusing!" % (len(set(pe_seen))))

    if pe_seen.count(1) != pe_seen.count(2):
        raise ValueError("Uneven pairing of paired ends (are you missing a file)? PE1 count: %d, PE2 count: %d" %
                         (pe_seen.count(1), pe_seen.count(2)))

    find_pe_logger.info("Files seen: %d" % (len(all_files)))
    find_pe_logger.info("Samples seen: %d" % (len(set(prefixes_seen))))

    merge_strategy = {}

    find_pe_logger.info("Sample IDs:")
    for prefix in sorted(set(prefixes_seen)):
        find_pe_logger.info("     %s" % (prefix))

    for file in file_dict.iterkeys():
        merge_strategy.setdefault(file_dict[file]['prefix'], []).append(file)

    if verbose:
        find_pe_logger.debug("Merge strategy is:")
        find_pe_logger.debug(pprint.pformat(merge_strategy))

    return merge_strategy

def run_bowtie2(paired_end_mapping, genome, output_path, disable_parallel=False):
    bowtie2_logger = _logshim.getLogger('run_bowtie2')

    # Import the config file to get genome locations
    config = _script_helpers.get_config()

    if disable_parallel:
        shell_job_runner = _script_helpers.ShellJobRunner(bowtie2_logger)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(bowtie2_logger, delay_seconds=60)

    for output_prefix, paired_ends in paired_end_mapping.iteritems():
        bowtie2_logger.info('Spawning niced process for bowtie2 on: %s' % (output_prefix))
        for filename in paired_ends:
            assert(" " not in filename)
            assert(";" not in filename)  # Vague sanity testing for input filenames
            bowtie2_logger.debug('    Input: %s' % (filename))

        # bowtie2 options:
        # --end-to-end: this is the default, but let's explicitly specify it
        # --sensitive: again, the default (consider switching to --fast?)
        # --no-unal: Suppress unaligned reads from the output .sam
        # --no-discordant: These are paired-end reads. We expect them to be non-discordant.
        # --mm: mmap MAP_SHARED (other processes can use our genome, cool!)
        # --met-stderr: Write metrics to stderr
        # --time: output the time things took
        # -x: target genome
        command = "bowtie2 --end-to-end --sensitive --no-unal --no-discordant --mm --met-stderr --time -x %s -1 %s -2 %s 2>%s | samtools view -bS - >%s"

        shell_job_runner.run(command % (config['bowtie2_genomes'][genome],
                                        paired_ends[0],
                                        paired_ends[1],
                                        output_path + "/" + output_prefix + ".bt2.log",
                                        output_path + "/" + output_prefix + ".bt2.bam"))

    shell_job_runner.finish()


def main():
    # Parse & interpret command line flags.
    parser = argparse.ArgumentParser(description='Given paired-end .fastq/.fastq.gz files, map to a genome.',
                                     epilog="Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at "
                                            "Washington University in St. Louis: https://gordonlab.wustl.edu.",
                                     usage='%(prog)s [options]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input-path', '-i', dest="input_path", metavar='input_dir/', type=str,
                        help='Input path.', required=True)
    parser.add_argument('--output-path', '-o', dest="output_path", metavar='output_dir/', type=str,
                        help='Output path.', required=True)
    parser.add_argument('--genome', '-g', dest="genome", metavar='genome', type=str,
                        choices=['mm9', 'mm10', 'hg18', 'hg19'], help='Genome to use for bowtie2', required=True)
    parser.add_argument('--no-parallel', '-np', dest="no_parallel", default=False, action='store_true',
                        help='Disable parallel job spawning.')


    parser.add_argument("--verbose", "-v", dest="verbose", default=False, action='store_true')

    parser.add_argument("--no-log", "-nl", dest="nolog", default=False, action='store_true',
                        help="Do not create a log file.")

    args = parser.parse_args()

    output_path = _script_helpers.setup_output_path(args.output_path)

    _logshim.startLogger(verbose=args.verbose, noFileLog=args.nolog, outPath=output_path)


    paired_end_mapping = find_paired_ends(args.input_path, verbose=args.verbose)

    run_bowtie2(paired_end_mapping, args.genome, output_path)


if __name__ == '__main__':
    main()
