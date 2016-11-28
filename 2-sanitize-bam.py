#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Do very simple cleaning of .bam files for ATAC-seq data pre-processing.
#
# This script:
# * Drops low mapping quality reads (<10)
# * Removes chrM mitochondrial reads
# * Removes PCR duplicates (using samtools)
# * Removes blacklisted regions (from ENCODE or custom blacklists)
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
import glob
import os
import tempfile

# A parameter needed by samtools to sort in-memory.
MAX_MEM = "50G"

# Load our config files
CONFIG = _script_helpers.get_config()


def large_filter_fixmate_and_sort(input_files, genome, output_path, disable_parallel=False):
    primary_logger = _logshim.getLogger('first_pass')

    output_suffix = ".tmp"

    if disable_parallel:  # Doesn't change parallelism in last samtools sort
        shell_job_runner = _script_helpers.ShellJobRunner(primary_logger)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(primary_logger, delay_seconds=60)

    # We do a few things here:
    #  - View only mapping quality >= 10
    #  - Remove chrM
    #  - Sort by name for fixmate
    #     We don't parallelize here (-@ #) because fixmate blocks & parallel seems to only help for compressed.
    #  - Fixmate (needed for rmdrup)
    #  - Resorted by position
    tempfiles = []
    for filename in input_files:
        primary_logger.debug('Working on: %s' % (filename))
        command = 'export LANG=C; %s view -h -q 10 %s | grep -vF "chrM" | %s view -u -b - | ' \
                  '%s sort -l 0 -n -m %s -T %s -O bam | %s fixmate -O bam - - | %s sort -@ 8 -m %s - %s'

        # A super evil user could modify TMPDIR and make this generate evil strings. That's evil.
        temporary_file = tempfile.mkstemp('.tmp.bam')
        tempfiles.append(temporary_file)

        shell_job_runner.run(command % (CONFIG['binaries']['samtools'],
                                        filename,
                                        CONFIG['binaries']['samtools'],
                                        CONFIG['binaries']['samtools'],
                                        MAX_MEM,
                                        temporary_file[1],
                                        CONFIG['binaries']['samtools'],
                                        CONFIG['binaries']['samtools'],
                                        MAX_MEM,
                                        output_path + "/" + os.path.basename(os.path.splitext(filename)[0]) + output_suffix))

    shell_job_runner.finish()

    # Clean up our temporary files.
    primary_logger.info('Removing temporary files ...')
    for fd, fname in tempfiles:
        os.close(fd)
        os.unlink(fname)

    primary_logger.info('First large stage complete! Saved as .tmp.bam for next stage.')


def rmdup_and_blacklist(input_files, genome, output_path, disable_parallel=False):
    primary_logger = _logshim.getLogger('rmdup_blacklist')

    if disable_parallel:
        shell_job_runner = _script_helpers.ShellJobRunner(primary_logger)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(primary_logger, delay_seconds=20)

    for filename in input_files:
        primary_logger.debug('Working on: %s' % (filename))
        # This is extremely fast and has minimal memory usage. Yay!
        # TODO: Allow adjustable windowing (-w %d) to blacklist larger/surrounding regions?
        command = "%s rmdup %s - 2>%s | %s window -abam - -b %s -v -w 0 > %s"

        shell_job_runner.run(command % (CONFIG['binaries']['samtools_legacy'],  # TODO: Update this when samtools is fixed.
                                        output_path + "/" + os.path.basename(os.path.splitext(filename)[0]) + '.tmp.bam',  # TODO: CLEAN THIS
                                        output_path + "/" + os.path.basename(os.path.splitext(filename)[0]) + '.srt.rmdup.bam.log',
                                        CONFIG['binaries']['bedtools'],
                                        os.path.dirname(os.path.realpath(__file__)) + '/' + CONFIG['blacklists'][genome],  # TODO: CLEAN THIS
                                        output_path + "/" + os.path.basename(os.path.splitext(filename)[0]) + '.srt.rmdup.bam'))

    shell_job_runner.finish()

    primary_logger.info('Removing temporary files from stage 1 ...')
    for filename in input_files:
        os.unlink(output_path + "/" + os.path.basename(os.path.splitext(filename)[0]) + '.tmp.bam')

    primary_logger.info('Completed rmdup and blacklist')


def main():
    # Parse & interpret command line flags.
    parser = argparse.ArgumentParser(description='Given input .bam files, fix matepairs, remove duplicates, blacklist bad'
                                                 'regions, and sort the output.',
                                     epilog="Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at "
                                            "Washington University in St. Louis: http://gordonlab.wustl.edu.",
                                     usage='%(prog)s [options]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input-path', '-i', dest="input_path", metavar='input_dir/', type=str,
                        help='Input path.', required=True)
    parser.add_argument('--output-path', '-o', dest="output_path", metavar='output_dir/', type=str,
                        help='Output path.', required=True)
    parser.add_argument('--genome', '-g', dest="genome", metavar='genome', type=str,
                        choices=['mm9', 'mm10', 'hg18', 'hg19'], help='Genome to use for blacklisting.', required=True)
    parser.add_argument('--no-parallel', '-np', dest="no_parallel", default=False, action='store_true',
                        help='Disable parallel job spawning.')

    parser.add_argument("--verbose", "-v", dest="verbose", default=False, action='store_true')

    parser.add_argument("--no-log", "-nl", dest="nolog", default=False, action='store_true',
                        help="Do not create a log file.")

    args = parser.parse_args()

    output_path = _script_helpers.setup_output_path(args.output_path)

    log_main = _logshim.startLogger(verbose=args.verbose, noFileLog=args.nolog, outPath=output_path)


    # Samtools requires a temp directory for sorting /sometimes/.
    # This seems to only matter if it exceeds the in-ram limits set by the MAX_MEM parameter.
    # Sanity check the /tmp directory has a bit of space.
    temppath = tempfile.gettempdir()
    s = os.statvfs(temppath)
    if ((s.f_bavail * s.f_frsize) / (1024 * 1024)) < 10000:  # ~10 G, not for any good reason though
        log_main.warn('Temp directory %s doesn\'t have a lot of free space!' % (temppath))


    input_files = glob.glob(args.input_path + "/*.bam")  # Take ALL of the .bams.

    large_filter_fixmate_and_sort(input_files, args.genome, output_path, disable_parallel=args.no_parallel)
    rmdup_and_blacklist(input_files, args.genome, output_path, disable_parallel=args.no_parallel)



if __name__ == '__main__':
    main()
