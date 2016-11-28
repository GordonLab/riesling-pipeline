#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script merges one sample across lanes.
#
# This script will *only* be useful if you have the same multiplexed sample loaded into multiple lanes of a flowcell.
#
# This concatenates the same index's pared-end files (L*_R*_* .fastq.gz) across multiple
# lanes into one set of PE files per-sample.
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


def fastq_map_predict(input_path, verbose=False):
    """
    Determine a sane .fastq muti-lane merge strategy.
    Fail if we can't merge correctly, if there are remaining files, etc.

    sample file name: Gordon-Ad2-11-AAGAGGCA-AAGAGGCA_S7_L001_R1_001.fastq.gz

    Args:
        input_path: An input path containing .fastq / .fastq.gz files
    Returns:
        A dict of mappings.
    """
    fastq_map_logger = _logshim.getLogger('fastq_map_predict')

    if not os.path.isdir(input_path):
        raise ValueError("Input must be a directory. You gave: %s" % (input_path))

    all_files = glob.glob(input_path + "/*_R*.fastq.gz")  # Ignore index files, must have _R in title
    all_files.extend(glob.glob(input_path + "/*_R*.fastq"))

    if len(all_files) == 0:
        raise ValueError("Input directory is empty!")

    # Given paired ends, we must always have an even number of input files.
    if len(all_files) % 2 != 0:
        raise ValueError("Input directory contains an odd number of files.")

    re_pattern = re.compile(r'^(.*)_L(\d+)_R(\d)_\d+(\.fastq|\.fastq\.gz)$')


    file_dict = OrderedDict()

    prefixes_seen = []
    lanes_seen = []
    pe_seen = []
    for file in sorted(all_files):
        if not os.access(file, os.R_OK):
            raise OSError("Cannot read file: %s" % (file))

        filename_only = file.rsplit('/', 1)[-1]
        result = re.match(re_pattern, filename_only)

        file_dict[file] = {'prefix': str(result.group(1)),
                           'L': int(result.group(2)),
                           'R': int(result.group(3))}

        prefixes_seen.append(file_dict[file]['prefix'])
        lanes_seen.append(file_dict[file]['L'])
        pe_seen.append(file_dict[file]['R'])


    # Sanity checking here. Missing files? Other oddities?
    if len(file_dict) % len(set(lanes_seen)) != 0:
        raise ValueError("Missing or extra file(s)? Saw %d lanes, and %d input files." %
                         (len(file_dict), len(set(lanes_seen))))

    if len(set(pe_seen)) != 2:
        raise ValueError("Saw %d paired ends, expecting exactly two. That's confusing!" % (len(set(pe_seen))))

    if pe_seen.count(1) != pe_seen.count(2):
        raise ValueError("Uneven pairing of paired ends (are you missing a file)? R1 count: %d, R2 count: %d" %
                         (pe_seen.count(1), pe_seen.count(2)))

    fastq_map_logger.info("Files seen: %d" % (len(all_files)))
    fastq_map_logger.info("Samples seen: %d" % (len(set(prefixes_seen))))
    fastq_map_logger.info("Lanes seen: %d" % (len(set(lanes_seen))))

    merge_strategy = {}

    fastq_map_logger.info("Sample IDs:")
    for prefix in sorted(set(prefixes_seen)):
        fastq_map_logger.info("     %s" % (prefix))

    for file in file_dict.iterkeys():
        merge_strategy.setdefault(file_dict[file]['prefix'] + ".PE" + str(file_dict[file]['R']), []).append(file)

    if verbose:
        fastq_map_logger.debug("Merge strategy is:")
        fastq_map_logger.debug(pprint.pformat(merge_strategy))

    return merge_strategy

def fastq_merge(merge_strategy, output_path, disable_parallel=False):
    """
    Concatenate multiple fastq files (from multiple lanes) into one.

    :param merge_strategy:
    :param output_path:
    :return:
    """
    merge_log = _logshim.getLogger('fastq_merge')

    if disable_parallel:
        shell_job_runner = _script_helpers.ShellJobRunner(merge_log)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(merge_log, delay_seconds=45)

    for merged_name, merge_inputs in merge_strategy.iteritems():
        merge_input_files = ' '.join(merge_inputs)
        merge_log.info('Spawning niced process to merge: %s' % (merged_name))
        for filename in merge_inputs:
            assert(" " not in filename)
            assert(";" not in filename)  # Vague sanity testing for input filenames
            merge_log.debug('    Input: %s' % (filename))

        # WARNING: Using shell has security implications! Don't work on untrusted input filenames.
        command = "zcat %s | gzip -1 > %s/%s.fastq.gz" % (merge_input_files, output_path, merged_name)

        shell_job_runner.run(command)

    shell_job_runner.finish()

    return True


def main():
    # Parse & interpret command line flags.
    parser = argparse.ArgumentParser(description='Intelligently merge fastq/fastq.gz files from an Illumina pipeline.'
                                     'Merges all L*_R*_* .fastq.gz files into one per sample.',
                                     epilog="Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at "
                                            "Washington University in St. Louis: http://gordonlab.wustl.edu.",
                                     usage='%(prog)s [options]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input-path', '-i', dest="input_path", metavar='input_dir/', type=str,
                        help='Input path.', required=True)
    parser.add_argument('--output-path', '-o', dest="output_path", metavar='output_dir/', type=str,
                        help='Output path.', required=True)
    parser.add_argument('--no-parallel', '-np', dest="no_parallel", default=False, action='store_true',
                        help='Disable parallel job spawning.')

    # parser.add_argument('--skip-stats', dest="skip_stats", action='store_true',
    #                    help='Skip statistics generation.', required=False)

    parser.add_argument("--verbose", "-v", dest="verbose", default=False, action='store_true')

    parser.add_argument("--no-log", "-nl", dest="nolog", default=False, action='store_true',
                        help="Do not create a log file.")

    args = parser.parse_args()

    output_path = _script_helpers.setup_output_path(args.output_path)

    _logshim.startLogger(verbose=args.verbose, noFileLog=args.nolog, outPath=output_path)

    # Our goal is to intelligently merge .fastq/.fastq.gz output from an Illumina run
    # The Illumina standard pipeline splits by barcode w/ semi-predictable filenames we can use, e.g.
    # IsoA-M1-CD4_S1_L001_I1_001.fastq.gz # index (discard)
    # IsoA-M1-CD4_S1_L001_R1_001.fastq.gz # end 1, lane 1
    # IsoA-M1-CD4_S1_L001_R2_001.fastq.gz # end 2, lane 2
    # IsoA-M1-CD4_S1_L002_I1_001.fastq.gz # index (discard), lane 2
    # IsoA-M1-CD4_S1_L002_R1_001.fastq.gz # end 1, lane 2
    # ...

    # TODO: Move some lower glob code up so we can test these functions
    merge_strategy = fastq_map_predict(args.input_path, verbose=args.verbose)

    fastq_merge(merge_strategy, args.output_path, disable_parallel=args.no_parallel)



if __name__ == '__main__':
    main()
