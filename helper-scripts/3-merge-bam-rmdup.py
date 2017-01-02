#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# In case you've sequenced the same sample multiple times, use this script to
# merge those samples together and then rmdup the pooled data.

# WARNING: This assumes you've used the pipeline up to this point, and have pre-sorted & fixmated .bams.
# Why use this script, when you could merge & rmdup by hand?
# * This just works (so you won't forget anything)
# * Automatic sane file naming
#
#
# Copyright (c) 2014-2016 Nick Semenkovich <semenko@alum.mit.edu>.
#   https://nick.semenkovich.com/
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
#   https://gordonlab.wustl.edu/
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

import _logshim
import _script_helpers
import argparse
import os

# Load our config files
CONFIG = _script_helpers.get_config()


def merge_and_rmdup(input_files, output_path, disable_parallel=False):
    primary_logger = _logshim.getLogger('merge_and_rmdup')

    # Sanity checks on the input files list
    assert(len(input_files) > 1)
    # Check files are readable
    for filename in input_files:
        if not os.access(filename, os.R_OK):
            primary_logger.fatal("Unable to read input files.")
            raise IOError

    output_file_name = '-AND-'.join([os.path.basename(os.path.splitext(filename)[0]) for filename in input_files])

    # Sanity check: maximum output filename length
    max_filename_length = os.statvfs(output_path).f_namemax
    if max_filename_length < 100:
        primary_logger.fatal("Cannot deal with short filename length limit. Maybe namemax is broken?")
        raise IOError

    if (len(output_file_name) + 10) > max_filename_length:  # roughly truncate filename for sanity.
        primary_logger.critical("Very long filename! Truncating!")
        output_file_name = output_file_name[:-20]  # Give us some extra room for downstream stuff?

    output_file_name += ".merged.bam"

    input_file_string = ' '.join(input_files)

    shell_job_runner = _script_helpers.ShellJobRunner(primary_logger)

    primary_logger.debug('Input file string: %s' % (input_file_string))
    primary_logger.debug('Working on merge as: %s' % (output_file_name))
    # This is pretty fast and has minimal memory usage. Yay!
    # We're probably re-rmduping some files if we're merging. That's ok since this is speedy.
    command = "%s merge -u - %s | %s rmdup - %s 2>%s"

    shell_job_runner.run(command % (CONFIG['binaries']['samtools'],
                                    input_file_string,
                                    CONFIG['binaries']['samtools_legacy'],  # TODO: Update this when samtools is fixed.
                                    output_path + "/" + output_file_name,
                                    output_path + "/" + os.path.basename(os.path.splitext(output_file_name)[0]) + '-rmdup.log'))

    shell_job_runner.finish()


    primary_logger.info('Merge and rmdup complete!')


def main():
    # Parse & interpret command line flags.
    parser = argparse.ArgumentParser(description='Pool multiple .bams together for the same sample.'
                                                 'Note: This is *only* necessary if you sequenced the same sample multiple times.',
                                     epilog="Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at "
                                            "Washington University in St. Louis: https://gordonlab.wustl.edu.",
                                     usage='%(prog)s [options]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input-files', '-i', dest="input_files", metavar='input_dir/', type=str,
                        help='Input files. (Not just a path!)', required=True, nargs='+')
    parser.add_argument('--output-path', '-o', dest="output_path", metavar='output_dir/', type=str,
                        help='Output path.', required=True)

    parser.add_argument("--verbose", "-v", dest="verbose", default=False, action='store_true')

    parser.add_argument("--no-log", "-nl", dest="nolog", default=False, action='store_true',
                        help="Do not create a log file.")

    args = parser.parse_args()

    output_path = _script_helpers.setup_output_path(args.output_path)

    _logshim.startLogger(verbose=args.verbose, noFileLog=args.nolog, outPath=output_path)

    merge_and_rmdup(args.input_files, output_path)



if __name__ == '__main__':
    main()
