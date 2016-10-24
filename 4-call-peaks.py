#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Nick Semenkovich <semenko@alum.mit.edu>.
#   https://nick.semenkovich.com/
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
#   http://gordonlab.wustl.edu/
#
# This software is released under the MIT License:
#  http://opensource.org/licenses/MIT
#
# Source: https://github.com/semenko/@@@@@@@@@@@@

# Compare multiple peak-calling strategies

# * PeakSeq
# * F-seq?
# * Hotspot??
# * GEM???


#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Nick Semenkovich <semenko@alum.mit.edu>.
#   https://nick.semenkovich.com/
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
#   http://gordonlab.wustl.edu/
#
# This software is released under the MIT License:
#  http://opensource.org/licenses/MIT
#
# Source: https://github.com/semenko/riesling-enhancer-prediction
#
# I gave good consideration to using a pipeline package (ruffus, etc.), but they
# seemed bloated or limited in their usefulness.


# Drop crappy very low (<@@@@@) from bowtie2
#   Remove ChrM reads >
#       Fix matepairs (samtools) >
#           Sort >
#               Remove duplicates >
#                   Remove blacklisted regions >
#                       Index >
#                           Generate Stats (% ChrM, etc.)

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


def generate_index(input_files, output_path, disable_parallel=False):
    """
    Many peak pickers want indexed .bams. Let's build indexes! (yay!)

    :param input_files:
    :param output_path:
    :param disable_parallel:
    :return:
    """
    primary_logger = _logshim.getLogger('index')

    if disable_parallel:
        shell_job_runner = _script_helpers.ShellJobRunner(primary_logger)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(primary_logger, delay_seconds=10)

    for filename in input_files:
        primary_logger.debug('Working on: %s' % (filename))
        command = "%s index %s"

        shell_job_runner.run(command % (CONFIG['binaries']['samtools'], filename))

    shell_job_runner.finish()



def run_macs14(input_files, output_path, genome, disable_parallel=False):
    macs14_log = _logshim.getLogger('run_macs14')

    macs14_log.info('Spawning MACS14 jobs...')

    if disable_parallel:
        shell_job_runner = _script_helpers.ShellJobRunner(macs14_log)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(macs14_log, delay_seconds=20)

    for filename in input_files:
        macs14_log.debug('Working on: %s' % (filename))

        # macs14 is old, but we try it anyway, since it's sometimes useful.
        # -t: input
        # -n: output name
        # -f: format
        # -g: genome
        # -p: pvalue for peak cutoff
        # --wig: save .wig outputs
        # --single-profile: make one single wiggle
        # --space=50: wiggle resolution (default: 10)
        #
        # Note: This CD hack is because MACS1.4 can't specify an output path :(
        command = "cd %s && %s -t %s -n %s -f BAM -g %s -p 1e-9 --wig --single-profile --space=50 2>%s"

        filename_without_extension = os.path.splitext(filename)[0] + '.macs14'

        shell_job_runner.run(command % (output_path,  # for cd hack
                                        'macs14',  # This must be pre-installed by the user. It's a big, complex package.
                                        os.getcwd() + '/' + filename,  # input file # TODO: Fix this path hack. MACS14 sucks and cannot specify an output path :/
                                        os.path.basename(filename_without_extension),
                                        genome,  # for genome size
                                        os.path.basename(filename_without_extension) + '.macs14.log'))

    shell_job_runner.finish()

    macs14_log.info('MACS14 peak calling complete.')


def run_macs2(input_files, output_path, genome, disable_parallel=False):
    macs2_log = _logshim.getLogger('run_macs2')

    macs2_log.info('Spawning MACS2 jobs...')

    if disable_parallel:
        shell_job_runner = _script_helpers.ShellJobRunner(macs2_log)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(macs2_log, delay_seconds=0.1)

    for filename in input_files:
        macs2_log.debug('Working on: %s' % (filename))

        # --bdg: generate .bed graph output
        # --nomodel: We'll be shifting manually!
        # --extsize 200: See long discussion at: @@@
        # --shift -100: As per above.
        # --slocal: Look at a local window of 20kb to build peak models
        # --keep-dup: We already removed duplicates with samtools.
        # TODO: Consider allowing tweaks to these settings with flags?
        command = "%s callpeak -t %s -n %s --outdir %s -g %s --bdg --nomodel --extsize 200 --shift -100 --slocal 20000 --llocal 20000 --keep-dup all 2>%s"

        filename_without_extension = os.path.splitext(filename)[0] + '.macs2'

        shell_job_runner.run(command % ('srun ~/.local/bin/macs2',  # This must be pre-installed by the user. It's a big, complex package.
                                        filename,  # input file
                                        os.path.basename(filename_without_extension),
                                        output_path,
                                        genome,  # for genome size, not sure this actually matters with nolambda/nomodel
                                        output_path + "/" + os.path.basename(filename_without_extension) + '.log'))

    shell_job_runner.finish()

    macs2_log.info('MACS2 peak calling complete.')


def run_homer(input_files, output_path, genome, disable_parallel=False):
    homer_log = _logshim.getLogger('run_homer')

    homer_log.info('Spawning Homer/findPeaks jobs...')

    raise NotImplementedError

    if disable_parallel:
        shell_job_runner = _script_helpers.ShellJobRunner(homer_log)
    else:
        shell_job_runner = _script_helpers.ShellJobRunner(homer_log, delay_seconds=20)

    for filename in input_files:
        homer_log.debug('Working on: %s' % (filename))

        # These are very sketchy, first-pass Homer settings.
        #
        command = ""

        filename_without_extension = os.path.splitext(filename)[0]

        shell_job_runner.run(command % (CONFIG['binaries']['findPeaks'],
                                        filename,  # input file
                                        filename_without_extension,
                                        output_path,
                                        genome,  # for genome size, not sure this actually matters with nolambda/nomodel
                                        output_path + "/" + filename_without_extension + '.log'))

    shell_job_runner.finish()

    homer_log.info('MACS2 peak calling complete.')

# findPeaks -i INPUT -o OUTPUT-FILE -size ??? -gsize SIZE -fragLength


def main():
    # Parse & interpret command line flags.
    parser = argparse.ArgumentParser(description='Run a number of standard peak calling algorithms for ATACseq data. '
                                                 'Expects de-duplicated, sorted, merged, ChrM-removed data.',
                                     epilog="Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at "
                                            "Washington University in St. Louis: http://gordonlab.wustl.edu.",
                                     usage='%(prog)s [options]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input-path', '-i', dest="input_path", metavar='input_dir/', type=str,
                        help='Input path (or a specific .bam file).', required=True)
    parser.add_argument('--output-path', '-o', dest="output_path", metavar='output_dir/', type=str,
                        help='Output path.', required=True)
    parser.add_argument('--genome', '-g', dest="genome", metavar='genome', type=str,
                        choices=['ms', 'mm', 'ce', 'dm'], help='Genome size to pass to MACS.', required=True)  # TODO: Consider using mm9/mm10, etc. for uniformity?
    parser.add_argument('--no-parallel', '-np', dest="no_parallel", default=False, action='store_true',
                        help='Disable parallel job spawning.')

    parser.add_argument('--skip-bam-indexing', dest="skip_bam_indexing", action='store_true',
                        help='Skip bam indexing (You must have generated indexes independently for peak callers to work!).', required=False)

    parser.add_argument('--skip-macs14', dest="skip_macs14", action='store_true',
                        help='Skip MACS v1.4 peak calling.', required=False)
    parser.add_argument('--skip-macs2', dest="skip_macs2", action='store_true',
                        help='Skip MACS v2 peak calling.', required=False)
    parser.add_argument('--skip-homer', dest="skip_homer", action='store_true',
                        help='Skip Homer (findPeaks) peak calling.', required=False)


    parser.add_argument("--verbose", "-v", dest="verbose", default=False, action='store_true')

    parser.add_argument("--no-log", "-nl", dest="nolog", default=False, action='store_true',
                        help="Do not create a log file.")

    args = parser.parse_args()

    output_path = _script_helpers.setup_output_path(args.output_path)

    log_main = _logshim.startLogger(verbose=args.verbose, noFileLog=args.nolog, outPath=output_path)

    input_files = _script_helpers.validate_input_files(args.input_path)


    # Generate BAM indexes
    if not args.skip_bam_indexing:
        generate_index(input_files, output_path, disable_parallel=args.no_parallel)
    else:
        log_main.warn("Skipping bam index .bai generation as requested.")
        log_main.warn("You must have generated these separately, otherwise peak callers will fail.")

    if not args.skip_macs14:
        # Start with old-school MACS 1.4
        run_macs14(input_files, output_path, args.genome, disable_parallel=args.no_parallel)

    if not args.skip_macs2:
        # Now new MACS 2
        # macs2 callpeak --nomodel -t $BAM -n $OUT --nolambda --keep-dup all --slocal 10000
        run_macs2(input_files, output_path, args.genome, disable_parallel=args.no_parallel)

    if not args.skip_homer and False:
        # NOT IMPLEMENTED (!)
        # Now the homer2 package
        run_homer(input_files, output_path, args.genome, disable_parallel=args.no_parallel)



if __name__ == '__main__':
    main()
