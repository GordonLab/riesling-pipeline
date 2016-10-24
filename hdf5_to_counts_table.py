#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script exists to convert .hdf5 files into counts tables readable by R.
#
# Unfortunately, we can't use the rhdf5 package, since it doesn't support the
# datatypes used by our .hdf5 files.
#
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

from __future__ import absolute_import, division, print_function, unicode_literals

__author__ = 'Nick Semenkovich <semenko@alum.mit.edu>'
__copyright__ = 'Gordon Lab at Washington University in St. Louis'
__license__ = 'MIT'
__version__ = '1.0.3'

import argparse
import csv
import fnmatch
import operator
import os
import pybedtools
import tables
import _logshim
import _script_helpers
from collections import deque, OrderedDict


CONFIG = _script_helpers.get_config()


# TODO: Modularize this function. This code is repeated in a *lot* of scripts.
def get_input_files(input_path):
    """
    Generate a list of all input files.

    :param input_files: A directory with .h5 files. (e.g. /tmp/)
    :return: a list of all .h5 files with absolute paths. (e.g. ['/tmp/a.h5'] )
    """
    if not os.path.isdir(input_path):
        raise ValueError("Input must be a directory. You gave: %s" % (input_path))

    # Adapted from:
    # https://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
    all_files = []
    for root, _, filenames in os.walk(input_path):
        for filename in fnmatch.filter(filenames, '*.h5'):
            all_files.append(os.path.join(root, filename))

    if len(all_files) == 0:
        raise ValueError("Input directory contains no .h5 files!")

    return all_files


def flatten_tsv(filename):
    """
    Flaten a TSV file -- parse and concatenate identical row names, by summing their values.
    """
    flatlog = _logshim.getLogger('flatten_tsv')

    flatlog.debug('Flattening input file: %s' % (filename))

    data_dict = OrderedDict()

    with open(filename, 'r') as tsv_ro_fh:
        tsv_input = csv.reader(tsv_ro_fh, delimiter=str("\t"))

        header = next(tsv_input, None)

        for row in tsv_input:
            row_key = row[0]
            these_row_values_as_int = map(int, row[1:])
            if row_key in data_dict:
                # Add the current row values to the existing values
                data_dict[row_key] = map(operator.add, data_dict[row_key], these_row_values_as_int)
            else:
                data_dict[row_key] = these_row_values_as_int

    # Write back the parsed dict
    with open(filename, 'wb') as tsv_rw_fh:
        tsv_writer = csv.writer(tsv_rw_fh, delimiter=str("\t"))
        tsv_writer.writerow(header)

        for key, val in data_dict.iteritems():
            tsv_writer.writerow([key] + val)



def parse_h5files(input_files, annotationBedTool, overwrite, flatten, density, normalized, sizescaled):
    h5logger = _logshim.getLogger('parse_h5files')

    assert(not (density and normalized))
    total_file_count = len(input_files)
    h5logger.info('Parsing a total of: %d file(s)' % (total_file_count))

    output_suffix_list = ['tsv']

    annotating_regions = False
    if annotationBedTool:
        annotating_regions = True
        output_suffix_list.append('annotated')

    if normalized:
        output_suffix_list.append('normalized')
    elif density:
        output_suffix_list.append('density')
    elif sizescaled:
        output_suffix_list.append('sizescaled')

    output_suffix = '.'.join(reversed(output_suffix_list))

    # Cache regions that we're annotating, maybe.
    region_annotation_cache = {}

    for this_file_count, file in enumerate(input_files):
        h5logger.info('\tParsing: %s (%d/%d)' % (file, this_file_count + 1, total_file_count))

        output_filename = file + '.' + output_suffix

        if not overwrite and os.path.isfile(output_filename):
            h5logger.warn('Skipping this .h5 as output .tsv already exists: %s' % (output_filename))
            continue

        # TODO: Modularize H5FD_CORE (the in-memory driver?)
        with tables.open_file(file, mode="r", driver="H5FD_CORE") as h5_object:
            assert(h5_object.title.startswith("bam liquidator genome read counts"))  # Some sanity checking
            assert(h5_object.root.file_names[0] == "*")

            bam_filename_header = h5_object.root.file_names[1:]
            bam_filename_header.insert(0, 'region')

            # Note: len(files) = len(file_names) - 1, since file_names has a 'wildcard' first entry.
            number_of_regions = int(len(h5_object.root.region_counts) / len(h5_object.root.files))

            # We expect this .h5 object's region_counts to contain:
            # /region_counts (Table(SIZE,)) 'region counts'
            #   description := {
            #   "file_key": UInt32Col(shape=(), dflt=0, pos=0),
            #   "chromosome": StringCol(itemsize=64, shape=(), dflt='', pos=1),
            #   "region_name": StringCol(itemsize=64, shape=(), dflt='', pos=2),
            #   "start": UInt64Col(shape=(), dflt=0, pos=3),
            #   "stop": UInt64Col(shape=(), dflt=0, pos=4),
            #   "strand": StringCol(itemsize=1, shape=(), dflt='', pos=5),
            #   "count": UInt64Col(shape=(), dflt=0, pos=6),
            #   "normalized_count": Float64Col(shape=(), dflt=0.0, pos=7)}
            #   byteorder := 'little'
            #   chunkshape := (NNN,)
            counts = h5_object.root.region_counts

            with open(output_filename, 'wb') as tsv_output:
                tsvwriter = csv.writer(tsv_output, delimiter=str("\t"))
                tsvwriter.writerow(bam_filename_header)

                if annotating_regions:
                    h5logger.debug('Generating .bed annotations from provided genome.')
                    region_to_gene = {}
                    # Perform one annotation rapidly for all regions in the .hdf5
                    hdf5_positions_only = []

                    for region_number in range(0, number_of_regions):
                        hdf5_positions_only.append(counts[region_number][1] + ' ' + str(counts[region_number][3]) + ' ' + str(counts[region_number][4]))

                    hdf5_positions_only_hashkey = ''.join(hdf5_positions_only)

                    if hdf5_positions_only_hashkey in region_annotation_cache:
                        # The genome doesn't change mid run, so we cache only on hdf5_positions
                        region_to_gene = region_annotation_cache[hdf5_positions_only_hashkey]
                        h5logger.debug('Annotation from cache.')
                    else:
                        hdf5_stub_bed = pybedtools.BedTool('\n'.join(hdf5_positions_only), from_string=True)

                        annotated_bed = hdf5_stub_bed.closest(annotationBedTool, t='first')

                        for locus in annotated_bed:
                            region_to_gene[locus.chrom + ':' + str(locus.start) + '-' + str(locus.end)] = locus.fields[11].split('"')[1]

                        region_annotation_cache[hdf5_positions_only_hashkey] = region_to_gene
                        h5logger.debug('Annotation completed.')


                # We're going to aggressively access the hdf5 at a bunch of fixed offsets.
                # rowarray = [counts[number_of_regions*0 + i], counts[number_of_regions*1 + i] + counts[number_of_regions*2 + i] ...]

                number_of_files = len(h5_object.root.files)
                working_deque = deque(maxlen=number_of_files + 1)

                # Here, we loop over every "region"/locus (every entry in the first column of the .tsv)
                # And then (within this loop) jump to each individual "file" (the hdf5 can contain multiple
                # separate samples) to build the data for every row.
                for region_number in range(0, number_of_regions):
                    # Prefix the row with chrN:bpSTART-pbEND e.g. chr4:100-2000
                    locus_name = counts[region_number][1] + ':' + str(counts[region_number][3]) + '-' + str(counts[region_number][4])

                    # Sanity checking, in case the input is nuts
                    feature_width = counts[region_number][4] - counts[region_number][3]
                    assert(feature_width > 0)

                    # DESeq2 requires each region have a unique name.
                    # You can either append a unique value, or aggregate identical loci.
                    # We address this later by re-opening and aggregating.
                    if annotating_regions:
                        working_deque.append(region_to_gene[locus_name])
                    else:
                        working_deque.append(locus_name)
                    #rowarray = [counts[region_number][1] + ':' + str(counts[region_number][3]) + '-' + str(counts[region_number][4])]

                    for file_number in range(0, number_of_files):
                        if normalized:
                            # Standard normalized (counts/mreads)
                            # bamliquidator gives us (counts/mreads)/width so we multiply by width
                            working_deque.append(int(counts[number_of_regions * file_number + region_number][7] * feature_width))
                        elif density:
                            # (counts/mreads)/width
                            # We upscale the fractional normalized count values by an arbitrary amount,
                            # because subsequent analyses like integers.
                            working_deque.append(int(counts[number_of_regions * file_number + region_number][7] * 10000))
                        elif sizescaled:
                            # counts/width
                            # We upscale the fractional normalized count values by an arbitrary amount,
                            # because subsequent analyses like integers.
                            working_deque.append(int(counts[number_of_regions * file_number + region_number][6] / feature_width * 100))
                        else:
                            working_deque.append(int(counts[number_of_regions * file_number + region_number][6]))

                    tsvwriter.writerow(working_deque)

            if flatten:
                flatten_tsv(output_filename)

    h5logger.info('Completed.')


def main():
    # Parse & interpret command line flags.
    parser = argparse.ArgumentParser(description='Convert hdf5 tables from bamliquidator format to CSV counts tables '
                                                 'for use in R and elsewhere. (Necessary as rhdf5 doesn\'t support our data structure.)',
                                     epilog="Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at "
                                            "Washington University in St. Louis: http://gordonlab.wustl.edu.",
                                     usage='%(prog)s [options]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input-path', '-i', dest="input_path", metavar='input_dir/', type=str,
                        help='Input path with .h5 files.',
                        required=True)

    parser.add_argument("--overwrite", dest="overwrite", default=False, action='store_true',
                        help='Regenerate and overwrite output .tsv files, even if they already exist.')

    parser.add_argument('--call-genes', dest="call_genes", default=False, action='store_true',
                        help='Instead of a .tsv (with positions as keys), make a .annotated.tsv with nearby genes.')

    parser.add_argument('--normalized', dest="normalized", default=False, action='store_true',
                        help='Store the normalized counts (counts/total reads) instead of the raw read counts.')

    parser.add_argument('--density', dest="density", default=False, action='store_true',
                        help='Store the width-normalized density (counts/total reads/region size) instead of the raw read counts.')

    parser.add_argument('--sizescaled', dest="sizescaled", default=False, action='store_true',
                        help='Store the size scaled counts (counts/feature size) instead of the raw read counts.')

    # Useful for EdgeR/DESeq2, etc. where every locus/position/gene-name must be unique.
    parser.add_argument('--flatten', dest="flatten", default=False, action='store_true',
                        help='Aggregate identical locus IDs and sum their values. '
                             'Think carefully before you sum non-normalized values!')


    genome_choices = sorted(CONFIG['gffs'].keys())
    parser.add_argument('--genome', '-g', dest="genome", metavar='genome', type=str, default=None,
                        choices=genome_choices, help='Genome to use for annotation, one of: %s' % (', '.join(genome_choices)), required=False)


    parser.add_argument("--verbose", "-v", dest="verbose", default=False, action='store_true')

    parser.add_argument("--no-log", "-nl", dest="nolog", default=False, action='store_true',
                        help="Do not create a log file.")

    args = parser.parse_args()

    if args.call_genes and not args.genome:
        parser.error('--genome is when requesting --call_genes')

    assert((args.density + args.normalized + args.sizescaled) <= 1)

    annotationBedTool = None
    if args.call_genes:
        genome_gff = CONFIG['gffs'][args.genome]
        assert(os.access(genome_gff, os.R_OK))
        annotationBedTool = pybedtools.BedTool(genome_gff)

    # Output path is input path. This also checks that the path is writeable.
    output_path = _script_helpers.setup_output_path(args.input_path)

    _logshim.startLogger(verbose=args.verbose, noFileLog=args.nolog, outPath=output_path)


    input_files = get_input_files(args.input_path)

    parse_h5files(input_files,
                  annotationBedTool=annotationBedTool,
                  overwrite=args.overwrite,
                  flatten=args.flatten,
                  density=args.density,
                  normalized=args.normalized,
                  sizescaled=args.sizescaled)



if __name__ == '__main__':
    main()
