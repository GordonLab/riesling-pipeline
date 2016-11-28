#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014-2016 Nick Semenkovich <semenko@alum.mit.edu>.
#   https://nick.semenkovich.com/
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
#   https://gordonlab.wustl.edu/
#
# This software is released under the MIT License:
#   http://opensource.org/licenses/MIT
#
# Source: https://github.com/GordonLab/riesling-pipeline

# Input:
#  - original BAM
#  - peaks from peak caller (.bed)
#  - bamliquidator output?
#
# Output:
#  - genes
#  - enh clusters
#  - stretch regions?
#  - evolutionary conservation at loci?
#  - some magic HTML output?
#  - filter really highly expressed / other identified stuff?

from __future__ import absolute_import, division, print_function, unicode_literals


__author__ = 'Nick Semenkovich <semenko@alum.mit.edu>'
__copyright__ = 'Gordon Lab at Washington University in St. Louis'
__license__ = 'MIT'
__version__ = '2.8.1'

import argparse
import copy
import cPickle
import _logshim
import csv
import os
import gzip
from operator import itemgetter
import pybedtools
from pybedtools import featurefuncs
import statistics
import tables
import _script_helpers
import tempfile
from collections import OrderedDict
from bamliquidatorbatch import bamliquidator_batch
from itertools import groupby

# Curve fitting / plotting of SE ops
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#import seaborn as sns

# TODO: Remove for prod
import IPython

# Load our config files
CONFIG = _script_helpers.get_config()


def mask_bed_file(input_bed, masking_bed):
    """
    Mask a bed file. Useful to remove hotspots, known bad regions, etc.

    :param input_bed: The input peak bed file (required!)
    :param masking_bed: An optional masking file (may be false, derived from args.mask_file)
    :return: A bedtool object of the peak file (masked or not, depending if masking_bed was a file)
    """
    mask_log = _logshim.getLogger('mask_bed_file')

    if masking_bed:
        mask_log.info('Masking input peaks .bed using: %s' % (masking_bed))
        mask_bed = pybedtools.BedTool(masking_bed)
        mask_log.info('  Mask feature count: %d' % (len(mask_bed)))
        original_peak_bedtool = pybedtools.BedTool(input_bed)
        mask_log.info('  Original peak count: %d' % (len(original_peak_bedtool)))
        peak_bedtool = original_peak_bedtool.subtract(mask_bed, A=True)  # A=Remove entire feature if ANY overlap.
    else:
        mask_log.info('No mask applied to input peaks .bed.')
        peak_bedtool = pybedtools.BedTool(input_bed)

    if len(peak_bedtool[0].name) == 0:
        raise InputError  # We need a name for each input .bed peak (make your bed four columns)

    # We'll remove the MACS_ prefix too.
    def _remove_macs_prefix(elt):
        if "macs" in str.lower(str(elt.name)):
            elt.name = elt.name.split("_")[-1]
        else:
            # There's no MACS prefix. Let's give each peak an ID corresponding to its row.
            _remove_macs_prefix.element_counter += 1
            elt.name = "%s" % (_remove_macs_prefix.element_counter)
        return elt
    _remove_macs_prefix.element_counter = 0

    # .each returns a generator, so we 'render' this with saveas()
    peak_bedtool_without_MACS_prefix = peak_bedtool.each(_remove_macs_prefix).saveas()

    # Make sure our peak names are all ints.
    #  This isn't strictly necessary. If you want to process non-int peak names,
    #  remove this and remove the clean_up_bed_feature_names function below.
    assert(int(elt.name) for elt in peak_bedtool_without_MACS_prefix) # Honestly, this just errors, but w/e.

    mask_log.info('  Peaks for analysis: %s' % len(peak_bedtool_without_MACS_prefix))
    assert(len(peak_bedtool_without_MACS_prefix) > 0)  # There's nothing left in the peak bed file! (Was it completely masked?)

    return peak_bedtool_without_MACS_prefix


def remove_transcription_start_sites(input_bedtool, ucsc_refseq_bedtool, ucsc_name_mapping, tss_window, remove_entire_feature=True):
    """
    Remove transcription start sites from a bedtool object.

    TSSes in ATAC have low SNR, generally clutter up what we want for super & stretch enhancers.

    (They may be useful for DESeq2-style analyses!)

    Used here to mask those peaks from MACS, etc.

    :param input_file:
    :param genome:
    :param tss_window: A double-sided window around which we exclude peaks.
    :return:
    """
    remove_tss_log = _logshim.getLogger('remove_TSSes')
    remove_tss_log.info('Removing TSSes and a window around them of: %d bp' % (tss_window))

    input_peak_count = len(input_bedtool)
    remove_tss_log.info('Number of input peaks: %d' % (input_peak_count))

    # TODO: Use name_mapping
    # The UCSC bedtool is the *entire* gene body. Let's scope that just to the TSS.
    ucsc_tsses_only = ucsc_refseq_bedtool.each(featurefuncs.TSS, upstream=tss_window, downstream=tss_window, add_to_name="_TSS")

    if remove_entire_feature:
        remove_tss_log.info('Removing *entire* peak if any overlap with TSS.')
    else:
        remove_tss_log.info('Only removing parts of peaks that overlap with TSS.')

    peak_bedtool_without_tsses = input_bedtool.subtract(ucsc_tsses_only, A=remove_entire_feature)  # A=remove *entire* feature if *any* overlap.

    remaining_peak_count = len(peak_bedtool_without_tsses)

    # This number may go UP if remove_entire_feature is true, since subtraction can then split features into two, etc.
    remove_tss_log.info('Approx. number of features removed: %d (%d remain)' % (abs(input_peak_count - remaining_peak_count),
                                                                                remaining_peak_count))

    # Fix the naming if we didn't remove whole features, since .subtract may cause peak name duplicates when features are split
    if not remove_entire_feature:
        peak_bedtool_without_tsses = fix_duplicate_bed_feature_names(peak_bedtool_without_tsses)

    assert(len(peak_bedtool_without_tsses) > 0)  # There's nothing left after removing TSSes! (Window too large?)

    return peak_bedtool_without_tsses


def load_refseq_ucsc_table(genome):
    """
    Load a .tsv RefSeq table into a BedTool object and name -> alternate_name dict.

    :param genome: The name of a genome in .config.yaml
    :return:
    """
    load_refseq_log = _logshim.getLogger('load_refseq')

    load_refseq_log.info('Loading RefSeq genes for: %s' % (genome))

    # TODO: Better error testing to assert genome in .config.yaml?

    name_mapping = {}
    # Generate a BedTool object from our RefSeq collection
    with gzip.open("%s/%s" % (_script_helpers.THISPATH, CONFIG['refseq'][genome]), 'rb') as tsvfile:
        bedStructure = []
        tsv_contents = csv.DictReader(tsvfile, dialect=csv.excel_tab)
        for item in tsv_contents:
            bedStructure.append((item['chrom'], item['txStart'], item['txEnd'], item['name'], 0, item['strand']))
            name_mapping[item['name']] = item['name2']
        ucsc_refseq_bedtool = pybedtools.BedTool(bedStructure)

    load_refseq_log.info('  Number of genes loaded: %d' % (len(name_mapping)))
    return (ucsc_refseq_bedtool, name_mapping)


def clean_up_bed_feature_names(peak_bedtool):
    """
    Convert bed feature names like:
    1_2_3_4_5 -> 1-5
     and
    1_2_3_5 -> 1-3_5

    Uses a clever approach, where we group index-value into bins. Assumes sorted input.
      Via https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list

    :param peak_bedtool: Any input bedtool with .name attributes
    :return:
    """
    clean_bed_log = _logshim.getLogger('clean_up_bed_feature_names')

    clean_bed_log.debug('Renaming bed features to make output sane.')

    def _clean_name(elt):
        name_ints = [int(x) for x in elt.name.split("_")]
        new_name = []
        for key, group in groupby(enumerate(name_ints), lambda (index, item): index - item):
            group = map(itemgetter(1), group)
            if len(group) > 1:
                new_name.append(str(group[0]) + '-' + str(group[-1]))
            else:
                new_name.append(str(group[0]))
        elt.name = '_'.join(new_name)
        return elt

    # .each returns a generator, so we 'render' this with saveas()
    cleaned_bedtool = peak_bedtool.each(_clean_name).saveas()

    return cleaned_bedtool


def fix_duplicate_bed_feature_names(peak_bedtool):
    """
    Due to TSS filtering via .subtract, we may end up with duplicate names.

    You could include this above, but I prefer the above function be generic.

    e.g.
    chr1 1-10 1 ->>
     chr1 1-5 1
     chr1 7-10 1

    Here, we rename these as:
    chr1 1-10 1 ->>
     chr1 1-5 1
     chr1 7-10 1_tssdupe1

    :param peak_bedtool:
    :return:
    """
    fix_dupe_logger = _logshim.getLogger('fix_duplicate_bed_feature_names')

    fix_dupe_logger.debug('Fixing duplicate names, since TSS removal can split peak ranges.')

    def _fix_dupes(elt):
        if str(elt.name) == _fix_dupes.last_name:
            _fix_dupes.this_dupe += 1
            elt.name += str('_tssdupe') + str(_fix_dupes.this_dupe)
        else:
            _fix_dupes.last_name = str(elt.name)
            _fix_dupes.this_dupe = 0
        return elt
    _fix_dupes.last_name = None
    _fix_dupes.this_dupe = 0

    # .each returns a generator, so we 'render' this with saveas()
    renamed_bedtool = peak_bedtool.each(_fix_dupes).saveas()

    return renamed_bedtool


def stitch_bed_using_window(peak_bedtool_without_tsses, output_path, stitch_override=None):
    """
    """
    stitch_bed_log = _logshim.getLogger('stitch_bed_using_window')

    stitch_bed_log.info('Stitching adjacent peaks in .bed...')

    if stitch_override is None:
        stitch_bed_log.info('Calculating dynamic stitching window.')
        stitch_window = calculate_stitching(peak_bedtool_without_tsses, output_path)
    else:
        stitch_window = stitch_override

    stitch_bed_log.info('Using a stitch window of: %d bp' % (stitch_window))
    stitched_peak_bedtool_without_tsses_unnamed = peak_bedtool_without_tsses.merge(d=stitch_window)

    #IPython.embed()

    # The .merge command dropped .name attributes. Let's add them back in, concatenating multiple names for merged peaks.
    #   c = use column 4 (.bed name)
    #   o = "collapse" together peak names
    #   delim = separate peak names with "_"

    # Is your script crashing here?
    # You have a broken version of BedTools and should update
    # See: https://github.com/arq5x/bedtools2/issues/175
    stitched_peak_bedtool_without_tsses = stitched_peak_bedtool_without_tsses_unnamed.map(peak_bedtool_without_tsses,
                                                                                          c=4, o='collapse', delim="_")

    # Let's clean up the names (1_2_3_4_7 -> 1-4_7)
    cleaned_stitched_peak_bedtool_without_tsses = clean_up_bed_feature_names(stitched_peak_bedtool_without_tsses)

    stitch_bed_log.info('  Remaining peaks after stitching: %d' % (len(cleaned_stitched_peak_bedtool_without_tsses)))

    return cleaned_stitched_peak_bedtool_without_tsses



def calculate_stitching(peak_bed_without_tsses, output_path):
    """
    This gets a little tricky. We're computing the ratio of window growth like this:

    |--A--|XX|--B--|   |--C--|
    [windowing by 2X]
    |----A----B----|   |--C--|

    Then, we compute the ratio:
        - gapped_peak_lengths mean (A + XX + B) + (C)
        - gapless_peak_lengths mean (A + B) + (C)

    That is, we're interested in the growth rate of the gaps (XX) to the cluster size.

    :param peak_bed_without_tsses:
    :param output_path:
    :return:
    """
    stitching_log = _logshim.getLogger('calculate_stitching')
    stitching_log.info('** Dynamically determining stitching window...')

    # Our maximum stitching window is 15KB, search in 500bp intervals
    merge_range = list(range(0, 15500, 500))

    cluster_counts = []
    gapped_peak_mean_lengths = []
    gapless_peak_mean_lengths = []

    for merging_window in merge_range:
        stitching_log.debug('  Clustering nearby regions with a window of: %d' % (merging_window))
        clustered_peaks = peak_bed_without_tsses.cluster(d=merging_window)

        # TODO: Clean up these list creations.
        # peaks_per_cluster = [0] * len(clustered_peaks)  # Zero fill a list.
        cluster_start_positions = [0] * len(clustered_peaks)
        cluster_gapless_lengths = []
        cluster_gapped_lengths = []

        cluster_id = None
        for peak in clustered_peaks:
            cluster_id = int(peak.fields[-1]) - 1  # The cluster number, assigned by bedtools. It starts at 1.
            # peaks_per_cluster[cluster_id] += 1
            # TODO: Also use .merge() from pybedtools
            if cluster_start_positions[cluster_id] == 0:
                cluster_start_positions[cluster_id] = peak.start
                cluster_gapless_lengths.append(0)
                cluster_gapped_lengths.append(0)
            cluster_gapless_lengths[cluster_id] = peak.end - cluster_start_positions[cluster_id]
            cluster_gapped_lengths[cluster_id] += len(peak)

        # Generate stats about this cluster generated by these merge parameters
        number_of_clusters = cluster_id + 1  # cluster_id is a zero-index counter.
        cluster_counts.append(number_of_clusters)

        gapless_peak_mean_lengths.append(statistics.mean(cluster_gapless_lengths))
        gapped_peak_mean_lengths.append(statistics.mean(cluster_gapped_lengths))

        # Consider using bedtools nuc to incorporate GC richness or other features?
        stitching_log.debug('  Unique peaks: %d' % (number_of_clusters))


    ### Select a (VERY ROUGH) minimum

    x_axis = cluster_counts
    y_axis = [x / y for x, y in zip(gapped_peak_mean_lengths, gapless_peak_mean_lengths)]
    # Compute differences between adjacent list values
    y_axis_derivative = [val - y_axis[i - 1] for i, val in enumerate(y_axis) if val != 1]  # 1 excludes our initial value

    # TODO: Use SciPy argrelextrema or something cleaner here?
    minimum_exclusion = 3  # Exclude the three first and three last stitching windows
    stitching_minimum_index, stitching_minimum_value = min(enumerate(y_axis_derivative[minimum_exclusion:-minimum_exclusion]),
                                                           key=itemgetter(1))
    stitch_window = merge_range[stitching_minimum_index + minimum_exclusion + 1]
    stitching_log.debug('Identified optimal stitching parameter: %d bp (d[gap]: %.4f)' % (stitch_window, stitching_minimum_value))


    ### Graph the output
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.title('Dynamic Windowing: Cutoff %s bp' % (stitch_window))
    plt.plot(x_axis, y_axis)
    plt.xlim([x_axis[0], x_axis[-1]])
    # Plot a vertical cutoff line
    plt.axvline(x=x_axis[merge_range.index(stitch_window)], color='red')
    plt.ylabel("Mean Peak-to-Peak Gap Size")

    # Plot the derivative
    plt.subplot(2, 1, 2)
    # We cut off one to get these to the same dimension
    try:
        plt.plot(x_axis[:-1], y_axis_derivative)
    except ValueError:
        # TODO: Investigate why this sometimes fails with edgecase inputs.
        stitching_log.warn('Dimensional error during graph generation. Graph abandoned.')
    plt.xlim([x_axis[0], x_axis[-1]])
    plt.axvline(x=x_axis[merge_range.index(stitch_window)], color='red')
    plt.xlabel("Number of Peaks")
    plt.ylabel("d(Mean Gap Size)")

    plt.savefig(output_path + '/dynstitch-window.png')
    plt.close()

    return stitch_window


def compute_read_density(stitched_peak_bedtool_without_tsses, input_bam, output_path, path_suffix=None, development_cache=False):
    """
    Use bamliquidator to compute the signal intensity (read density)

    :return:
    """
    density_log = _logshim.getLogger('compute_read_density')
    density_log.info('Computing read density in .bam using bamliquidator.')

    if path_suffix is not None:
        bamliquidator_output_path = _script_helpers.setup_output_path(output_path + '/bamliquidator' + path_suffix)
    else:
        bamliquidator_output_path = _script_helpers.setup_output_path(output_path + '/bamliquidator/')


    if development_cache:
        try:
            with open(bamliquidator_output_path + '/pickle.cache', 'rb') as cache_fh:
                cached_data = cPickle.load(cache_fh)
                density_log.warn('Using cache! ONLY use this for development purposes!')
            return cached_data
        except (IOError, cPickle.UnpicklingError):
            density_log.warn('  Cache corrupt or doesn\'t yet exist. Ignoring.')


    bamliquidator_extension_window = 200
    density_log.info('Running bamliquidator_batch using an extension window of: %d bp' % (bamliquidator_extension_window))
    liquidator = bamliquidator_batch.RegionLiquidator(stitched_peak_bedtool_without_tsses.fn,
                                                      bamliquidator_output_path, input_bam, region_format='bed',
                                                      extension=bamliquidator_extension_window, sense='.')

    # TODO: Put this behind a flag to create debug matrix.gff, or if people are interested in output.
    if False:
        bamliquidator_batch.write_bamToGff_matrix(os.path.join(bamliquidator_output_path, "matrix.gff"), liquidator.counts_file_path)


    density_log.info('Loading and sanity-checking bamliquidator data.')
    # We need to map the bamliquidator output (density at a region) back to genes / nearby features.

    bamliquidator_data = OrderedDict()
    with tables.openFile(bamliquidator_output_path + "/counts.h5", driver="H5FD_CORE") as bamliquidator_output_table:  # H5FD_CORE = in-memory driver
        assert(bamliquidator_output_table.root.files.nrows == 1)  # Check we only have one .bam input represented in the .h5

        valid_chromosomes = set(['chrX'] + ['chr' + str(i) for i in range(1, 30)])  # Valid chromosomes (not ChrY, not ChrM)

        skipped_due_to_invalid_regions = 0
        for row in bamliquidator_output_table.root.region_counts:
            # description := {
            # "file_key": UInt32Col(shape=(), dflt=0, pos=0),
            # "chromosome": StringCol(itemsize=16, shape=(), dflt='', pos=1),
            # "region_name": StringCol(itemsize=64, shape=(), dflt='', pos=2),
            # "start": UInt64Col(shape=(), dflt=0, pos=3),
            # "stop": UInt64Col(shape=(), dflt=0, pos=4),
            # "strand": StringCol(itemsize=1, shape=(), dflt='', pos=5),
            # "count": UInt64Col(shape=(), dflt=0, pos=6),
            # "normalized_count": Float64Col(shape=(), dflt=0.0, pos=7)}
            # byteorder := 'little'
            # chunkshape := (560,)
            assert(row['file_key'] == 1)
            if row['chromosome'] in valid_chromosomes:
                assert(row['region_name'] not in bamliquidator_data)  # Peak IDs should've been unique in the original .bed. Otherwise we'll collide here.
                assert(row['stop'] > row['start'])
                bamliquidator_data[row['region_name']] = {'chromosome': row['chromosome'], 'strand': row['strand'],
                                                          'start': row['start'], 'stop': row['stop'],
                                                          'count': row['count'], 'normalized_count': row['normalized_count'],
                                                          'width_normalized_count': row['normalized_count'] * (row['stop'] - row['start'])}
            else:
                skipped_due_to_invalid_regions += 1

        density_log.info('Bamliquidator entries skipped due to invalid chromosome: %d (%.2f%%)' %
                         (skipped_due_to_invalid_regions,
                          (float(skipped_due_to_invalid_regions) / bamliquidator_output_table.root.region_counts.nrows) * 100))

        if skipped_due_to_invalid_regions > 50:
            density_log.warn('High number of skipped chromosomes! Did you filter the input data?')

    if development_cache:
        density_log.warn('Writing to development cache! Do not use in production!')
        with open(bamliquidator_output_path + '/pickle.cache', 'wb') as cache_fh:
            cPickle.dump(bamliquidator_data, cache_fh, cPickle.HIGHEST_PROTOCOL)

    return bamliquidator_data


def get_stats(any_list):
    """
    Return basic stats about a list
    :param any_list: a list
    :return: mean, median, stdev, variance
    """
    mean = statistics.mean(any_list)
    median = statistics.median(any_list)
    stddev = statistics.stdev(any_list)
    variance = statistics.variance(any_list)

    return mean, median, stddev, variance


def call_enhancers(bamliquidator_data, output_path, name_prefix):
    """
    This a top 10% enhancer population. Use the R code for the Young-style tangent line cutoff.

    :param bamliquidator_data:
    :param output_prefix:
    :param name_prefix:
    :return:
    """
    super_log = _logshim.getLogger('call_enhancers')
    super_log.info('** Predicting enhancers...')
    super_log.info('  Note: This uses 10% signal bounds, rather than the tangent cutoff. Use .R code for tangent.')

    raw_width_normalized_counts = [item['width_normalized_count'] for item in bamliquidator_data.values()]
    assert(min(raw_width_normalized_counts) >= 0)  # Only in the control-adjusted data, below.

    mean, median, stddev, variance = get_stats(raw_width_normalized_counts)

    num_enhancers = len(raw_width_normalized_counts)

    super_log.info('Total enhancer loci: %d' % (num_enhancers))
    super_log.info('  Mean: %.2f, Median: %.2f, Stdev: %.2f, Variance: %.2f' % (mean, median, stddev, variance))

    # Remember Chebyshev? Because man, I sure don't.
    # 1.281 = 80% w/i CI, 20% outside the CI (it's double sided, so this selects the top 10%)
    # 2.575 = 99% w/i CI -> 1% outside -> top 0.5% (!)
    super_enhancer_cutoff = mean + (1.281 * stddev)

    super_log.warn('Curve fitting not working! Using SD bounds.')
    super_log.info('  Choosing enhancer cutoff of: %.2f' % (super_enhancer_cutoff))


    # We know SEs can make up enormous amounts of input signal (sometimes >40%).
    # Let's use the Young group suggestion of: plot & approximate tangent
    # So we'll do this by:
    # * Plot all points, ranked by signal
    # * Fit an exponential
    # * Find the tangent with slope 1, D(fit curve) = 1
    # TODO: Implement R code here? Or leave in R?

    # bamliquidator_data keys, descending sort based on width_normalized_count
    sorted_bamliquidator_keys = sorted(bamliquidator_data, key=lambda k: bamliquidator_data[k]['width_normalized_count'])

    x = range(0, num_enhancers)
    y_width_normalized_counts = [bamliquidator_data[key]['width_normalized_count'] for key in sorted_bamliquidator_keys]

    # y_raw_counts = [bamliquidator_data[key]['count'] for key in sorted_bamliquidator_keys]

    # When does the count go over our cutoff?
    cutoff_super_vs_traditional = None
    for i, value in enumerate(y_width_normalized_counts):
        if value > super_enhancer_cutoff:
            cutoff_super_vs_traditional = i
            break
    assert(cutoff_super_vs_traditional is not None)
    super_log.debug('SE vs TE cutoff at %d' % (cutoff_super_vs_traditional))

    largest_val = max([abs(elt) for elt in y_width_normalized_counts])
    y_rescaled_1 = [(float(i) / largest_val) * 1 for i in y_width_normalized_counts]

    plt.figure()
    plt.title('Enhancer Signal Distribution')
    plt.xlabel('Enhancers Ranked by Signal')
    plt.ylabel('Normalized Signal')
    plt.plot(y_rescaled_1, label='Enhancer Signal', color='black')

    plt.xlim(0, num_enhancers)
    plt.ylim(0, 1)

    # Plot a vertical cutoff line
    plt.axvline(x=cutoff_super_vs_traditional, color='green')

    # Shade the TE area
    plt.fill_between(x[:cutoff_super_vs_traditional], y_rescaled_1[:cutoff_super_vs_traditional], facecolor='blue', alpha=0.5)

    # Shade the SE area
    plt.fill_between(x[cutoff_super_vs_traditional:], y_rescaled_1[cutoff_super_vs_traditional:], facecolor='red', alpha=0.5)

    ### Stats for the .text()
    total_normalized_signal = sum(y_width_normalized_counts)
    se_signal = sum(y_width_normalized_counts[cutoff_super_vs_traditional:])
    te_signal = sum(y_width_normalized_counts[:cutoff_super_vs_traditional])


    # Add some figure text
    # Our xpos is 10% of # enhancers to get a bit offset from the axis
    plt.text(num_enhancers * .1, 0.6,
             'Agg. Signal: %d\n\n#SE: %d (%.2f%%)\n  Signal: %d (%.2f%%)\n\n#TE: %d (%.2f%%)\n  Signal: %d (%.2f%%)' %
             (total_normalized_signal,
              num_enhancers - cutoff_super_vs_traditional,
              (1 - (float(cutoff_super_vs_traditional) / num_enhancers)) * 100,
              se_signal,
              float(se_signal) / total_normalized_signal * 100,

              cutoff_super_vs_traditional,
              (float(cutoff_super_vs_traditional) / num_enhancers) * 100,
              te_signal,
              float(te_signal) / total_normalized_signal * 100),
             style='italic', bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})

    plt.legend(loc=2)  # Legend in 'upper left'
    plt.savefig(output_path + '/' + name_prefix + '-super-vs-traditional-enhancers.png')
    plt.close()

    se_peak_ids = sorted_bamliquidator_keys[cutoff_super_vs_traditional:]
    te_peak_ids = sorted_bamliquidator_keys[:cutoff_super_vs_traditional]

    return(se_peak_ids, te_peak_ids)



def call_enhancers_with_control(bamliquidator_data, output_path, name_prefix):
    """
    Same as the above "call_enhancers" except this assumes a sigmoidal-esque function.

    :param bamliquidator_data:
    :param sorted_bamliquidator_keys:
    :param output_path:
    :return:
    """
    super_log = _logshim.getLogger('call_enhancers_with_control')
    super_log.info('** Predicting enhancers with control...')

    raw_width_normalized_counts = [item['width_normalized_count'] for item in bamliquidator_data.values()]
    assert(min(raw_width_normalized_counts) >= 0)  # Simple sanity checking.

    mean, median, stddev, variance = get_stats(raw_width_normalized_counts)

    num_enhancers = len(raw_width_normalized_counts)

    super_log.info('Total enhancer loci: %d' % (num_enhancers))
    super_log.info('  Mean: %.2f, Median: %.2f, Stdev: %.2f, Variance: %.2f' % (mean, median, stddev, variance))

    # Remember Chebyshev? Because man, I sure don't.
    # 1.281 = 80% w/i CI, 20% outside the CI (it's double sided, so this selects the top 10%)
    # 2.575 = 99% w/i CI -> 1% outside -> top 0.5% (!)
    super_enhancer_cutoff = mean + (1.281 * stddev)
    neg_super_enhancer_cutoff = mean - (1.281 * stddev)

    super_log.warn('Curve fitting not working! Using SD bounds.')
    super_log.info('  Global enhancer cutoff of: %.2f' % (super_enhancer_cutoff))
    super_log.info('  Global neg. enhancer cutoff of: %.2f' % (neg_super_enhancer_cutoff))

    # We know SEs can make up enormous amounts of input signal (sometimes >40%).
    # Let's use the Young group suggestion of: plot & approximate tangent
    # So we'll do this by:
    # * Plot all points, ranked by signal
    # * Fit an exponential
    # * Find the tangent with slope 1, D(fit curve) = 1
    # TODO: Implement R code here? Or leave in R?

    # bamliquidator_data keys, descending sort based on normalized_count
    sorted_bamliquidator_keys = sorted(bamliquidator_data, key=lambda k: bamliquidator_data[k]['width_normalized_count'])

    x = range(0, num_enhancers)
    y_width_normalized_counts = [bamliquidator_data[key]['width_normalized_count'] for key in sorted_bamliquidator_keys]
    # The index of y_normalized_counts where we switch from negative (adjusted) to zero/positive.

    # TODO: FIX THIS. No longer neg/pos after width scaling.
    zero_point = min(enumerate(y_width_normalized_counts), key=lambda i: abs(i[1] - 1))[0]  # both sides of 1



    ### For the control adjusted, we determine cutoffs for both the + and - sides independently.
    pos_mean, pos_median, pos_stddev, pos_variance = get_stats(y_width_normalized_counts[zero_point:])
    neg_mean, neg_median, neg_stddev, neg_variance = get_stats(y_width_normalized_counts[:zero_point])

    # Override the above variables
    super_enhancer_cutoff = pos_mean + (1.281 * pos_stddev)
    neg_super_enhancer_cutoff = neg_mean - (1.281 * neg_stddev)
    super_log.info('+ Mean: %.2f, Median: %.2f, Stdev: %.2f, Variance: %.2f' % (pos_mean, pos_median, pos_stddev, pos_variance))
    super_log.info('+  POS enhancer cutoff of: %.2f' % (super_enhancer_cutoff))
    super_log.info('- Mean: %.2f, Median: %.2f, Stdev: %.2f, Variance: %.2f' % (neg_mean, neg_median, neg_stddev, neg_variance))
    super_log.info('-  NEG enhancer cutoff of: %.2f' % (neg_super_enhancer_cutoff))


    # y_raw_counts = [bamliquidator_data[key]['count'] for key in sorted_bamliquidator_keys]

    # When does the count go over our cutoff?
    cutoff_super_vs_traditional = None
    for i, value in enumerate(y_width_normalized_counts):
        if value > super_enhancer_cutoff:
            cutoff_super_vs_traditional = i
            break
    assert(cutoff_super_vs_traditional is not None)
    super_log.debug('SE vs TE cutoff at %d' % (cutoff_super_vs_traditional))

    neg_cutoff_super_vs_traditional = None
    for i, value in enumerate(y_width_normalized_counts):
        if value > neg_super_enhancer_cutoff:
            neg_cutoff_super_vs_traditional = i
            break
    assert(neg_cutoff_super_vs_traditional is not None)
    super_log.debug('Negative SE vs TE cutoff at %d' % (neg_cutoff_super_vs_traditional))

    largest_val = max([abs(elt) for elt in y_width_normalized_counts])
    y_rescaled_1 = [(float(i) / largest_val) * 1 for i in y_width_normalized_counts]

    ##### TESTING
    pos_largest_val = max([abs(elt) for elt in y_width_normalized_counts[zero_point:]])
    neg_largest_val = max([abs(elt) for elt in y_width_normalized_counts[:zero_point]])

    pos_y_rescaled_1 = [(float(i) / pos_largest_val) * 1 for i in y_width_normalized_counts[zero_point:]]
    neg_y_rescaled_1 = [(float(i) / neg_largest_val) * 1 for i in y_width_normalized_counts[:zero_point]]

    y_rescaled_1 = neg_y_rescaled_1 + pos_y_rescaled_1
    # IPython.embed()
    ##### TESTING

    plt.figure()
    plt.title('Enhancer Signal (Control Adjusted, Agg. Signal: %d)' % (sum(y_width_normalized_counts)))
    plt.xlabel('Enhancers Ranked by Signal')
    plt.ylabel('Normalized Signal')
    plt.plot(y_rescaled_1, color='black')  # label='Enhancer Signal'

    plt.xlim(0, num_enhancers)
    plt.ylim(-1, 1)

    # Plot a vertical line for zero point
    plt.axvline(x=zero_point, linewidth=1.0, dashes=[8, 4], color='black')  # 0.9 = gray

    # Plot a vertical cutoff line for + SEs
    plt.axvline(x=cutoff_super_vs_traditional, linewidth=1.0, color='green')

    # Plot a vertical cutoff line for - SEs
    plt.axvline(x=neg_cutoff_super_vs_traditional, linewidth=1.0, color='green')

    # Shade the + TE area
    plt.fill_between(x[zero_point:cutoff_super_vs_traditional], y_rescaled_1[zero_point:cutoff_super_vs_traditional], facecolor='blue', alpha=0.5)
    # Shade the + SE area
    plt.fill_between(x[cutoff_super_vs_traditional:], y_rescaled_1[cutoff_super_vs_traditional:], facecolor='red', alpha=0.5)
    # Shade the - TE area
    plt.fill_between(x[neg_cutoff_super_vs_traditional:zero_point], y_rescaled_1[neg_cutoff_super_vs_traditional:zero_point], facecolor='blue', alpha=0.5)
    # Shade the - SE area
    plt.fill_between(x[:neg_cutoff_super_vs_traditional], y_rescaled_1[:neg_cutoff_super_vs_traditional], facecolor='red', alpha=0.5)


    ## Some stats for the .text()
    abs_normalized_counts = map(abs, y_width_normalized_counts)
    total_abs_normalized_counts = sum(abs_normalized_counts)

    se_signal = sum(abs_normalized_counts[cutoff_super_vs_traditional:])
    te_signal = sum(abs_normalized_counts[zero_point:cutoff_super_vs_traditional])

    neg_se_signal = sum(abs_normalized_counts[:neg_cutoff_super_vs_traditional])
    neg_te_signal = sum(abs_normalized_counts[neg_cutoff_super_vs_traditional:zero_point])

    # pos_normalized_counts = sum(abs_normalized_counts[zero_point:])
    # neg_normalized_counts = sum(abs_normalized_counts[:zero_point])

    # Add some figure text
    # Our xpos is 15% of # enhancers to get a bit offset from the axis
    plt.text(num_enhancers * .15, 0.2,
             '#PosSE: %d (%.2f%%)\n  Signal: %d (%.2f%%)\n\n#PosTE: %d (%.2f%%)\n  Signal: %d (%.2f%%)' %
             (num_enhancers - cutoff_super_vs_traditional,
              (1 - (float(cutoff_super_vs_traditional) / num_enhancers)) * 100,
              se_signal,
              float(se_signal) / total_abs_normalized_counts * 100,

              cutoff_super_vs_traditional - zero_point,
              (float(cutoff_super_vs_traditional - zero_point) / num_enhancers) * 100,
              te_signal,
              float(te_signal) / total_abs_normalized_counts * 100),
             bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})  # style='italic',

    plt.text(num_enhancers * .15, -0.8,
             '#NegSE: %d (%.2f%%)\n  Signal: %d (%.2f%%)\n\n#NegTE: %d (%.2f%%)\n  Signal: %d (%.2f%%)' %
             (neg_cutoff_super_vs_traditional,
              (float(neg_cutoff_super_vs_traditional) / num_enhancers) * 100,
              neg_se_signal,
              float(neg_se_signal) / total_abs_normalized_counts * 100,

              zero_point - neg_cutoff_super_vs_traditional,
              (float(zero_point - neg_cutoff_super_vs_traditional) / num_enhancers) * 100,
              neg_te_signal,
              float(neg_te_signal) / total_abs_normalized_counts * 100),
             bbox={'facecolor': 'red', 'alpha': 0.1, 'pad': 10})  # style='italic',

    plt.legend(loc=2)  # Legend in 'upper left'
    plt.savefig(output_path + '/' + name_prefix + '-super-vs-traditional-enhancers.png')
    plt.close()

    te_peak_ids = sorted_bamliquidator_keys[zero_point:cutoff_super_vs_traditional]
    se_peak_ids = sorted_bamliquidator_keys[cutoff_super_vs_traditional:]

    neg_te_peak_ids = sorted_bamliquidator_keys[neg_cutoff_super_vs_traditional:zero_point]
    neg_se_peak_ids = sorted_bamliquidator_keys[:neg_cutoff_super_vs_traditional]

    return(se_peak_ids, te_peak_ids, neg_se_peak_ids, neg_te_peak_ids)


def write_riesling_data_to_tsv(bamliquidator_data, output_file=None, keys=None):
    if output_file is None:
        raise NotImplementedError

    # TODO: Consider rearchitecting this around OrderedDict and DictWriter?
    if keys is not None:
        assert(isinstance(keys, list))

        with open(output_file, 'wb') as tsvfile:
            writer = csv.writer(tsvfile, dialect='excel-tab')
            writer.writerow(['peak_id', 'start', 'stop', 'width_normalized_count', 'normalized_count'])
            for key in keys:
                data = bamliquidator_data[key]
                writer.writerow([key, data['start'], data['stop'], data['chromosome'], data['width_normalized_count'], data['normalized_count']])
    else:
        raise NotImplementedError


def peak_id_to_bed(bamliquidator_data, peak_ids=None):
    """
    Generate a bed object from the bamliquidator data structure.

    ** Assumes bamliquidator_data is sorted. Does not do any sorting itself! **

    Very useful for then exporting files, doing usual .bed intersections, etc.

    :param bamliquidator_data:
    :param peak_ids:
    :return:
    """
    bed_creator_log = _logshim.getLogger('peak_id_to_bed')

    bedStructure = []

    if peak_ids is None:
        bed_creator_log.info('Converting entire bamliquidator structure to .bed object.')
        # No peaks provided? Dump the whole data structure.
        # TODO: OrderedDict the bamliqudiator_data or otherwise sort? ***
        for peak_id, peak_data in bamliquidator_data.iteritems():
            bedStructure.append((peak_data['chromosome'], peak_data['start'], peak_data['stop'], peak_id, peak_data['width_normalized_count'], peak_data['strand'], peak_data['normalized_count'], peak_data['stop'] - peak_data['start']))
    else:
        # Only include the selected peak_ids in a .bed object
        assert(isinstance(peak_ids, list))
        bed_creator_log.info('Creating .bed object from %d peaks in bamliquidator data.' % (len(peak_ids)))

        for peak_id in peak_ids:
            peak_data = bamliquidator_data[peak_id]
            bedStructure.append((peak_data['chromosome'], peak_data['start'], peak_data['stop'], peak_id, peak_data['width_normalized_count'], peak_data['strand'], peak_data['normalized_count'], peak_data['stop'] - peak_data['start']))

    bedtool_object = pybedtools.BedTool(bedStructure)
    return bedtool_object


def multipeak_stretch_to_bed(bamliquidator_data, filtered_stretch_enhancer_clusters):
    """
    Generate a bed object from the bamliquidator data structure.

    This works on the *filtered stretch enhancer cluster data*


    :param bamliquidator_data:
    :param filtered_stretch_enhancer_clusters: a list of [ ('8784', '8785', '8786'): 18.33837464700345, .... ]
    :return:
    """
    stitched_bed_creator_log = _logshim.getLogger('peak_id_to_bed')

    stitched_bed_creator_log.debug('Creating .bed object from %d stretches.' % (len(filtered_stretch_enhancer_clusters)))

    bedStructure = []
    for peak_cluster, cluster_score in filtered_stretch_enhancer_clusters.iteritems():

        first_peak_id = peak_cluster[0]
        last_peak_id = peak_cluster[-1]

        merged_peak = {
            'chromosome': bamliquidator_data[first_peak_id]['chromosome'],
            'start': bamliquidator_data[first_peak_id]['start'],
            'stop': bamliquidator_data[last_peak_id]['stop'],
            'strand': bamliquidator_data[first_peak_id]['strand']
        }

        aggregated_count = 0
        for peak_id in peak_cluster:
            # Lots of sanity checks
            assert(bamliquidator_data[peak_id]['chromosome'] == merged_peak['chromosome'])
            assert(bamliquidator_data[peak_id]['start'] >= merged_peak['start'])
            assert(bamliquidator_data[peak_id]['stop'] <= merged_peak['stop'])
            assert(bamliquidator_data[peak_id]['strand'] == merged_peak['strand'])
            aggregated_count += bamliquidator_data[peak_id]['normalized_count']  # TODO: Use raw count?

        pseudo_peak_id = '_'.join(peak_cluster)
        bedStructure.append((merged_peak['chromosome'], merged_peak['start'], merged_peak['stop'], pseudo_peak_id, aggregated_count, merged_peak['strand']))

    bedtool_object = pybedtools.BedTool(bedStructure)
    return bedtool_object


def update_signal_using_control(sample_bamliquidator_data, control_bamliquidator_data):
    """
    Remove background (control/comparison group) signal from SE/TEs.
    """
    update_signal_log = _logshim.getLogger('update_signal_using_control')

    update_signal_log.info('Updating raw and normalized counts using control.')
    assert(len(sample_bamliquidator_data) == len(control_bamliquidator_data))

    updated_dict = copy.deepcopy(sample_bamliquidator_data)

    for control_id, control_data in control_bamliquidator_data.iteritems():
        # WARNING: You can't zip() in sample/control data: dict's aren't guaranteed to be ordered!

        updated_dict[control_id]['width_normalized_count'] /= control_data['width_normalized_count'] + 1
#        updated_dict[control_id]['width_normalized_count'] *= (control_data['stop'] - control_data['start'])

        updated_dict[control_id]['normalized_count'] /= control_data['normalized_count'] + .0000000001  # TODO: Should we subtract here? Other options to adust? Recalculate?

        ################################ WILD TESTING CODE HERE
        # Trying some wild division instead
        #ncount_1 = updated_dict[control_id]['normalized_count']
        #ncount_2 = control_data['normalized_count']
        #try:
        #    if ncount_1 >= ncount_2:
        #        updated_dict[control_id]['normalized_count'] = (ncount_1/(ncount_2+0.000000001)) - 1
        #    else:
        #        updated_dict[control_id]['normalized_count'] = -((ncount_2/(ncount_1+0.000000001)) - 1)
        #except ZeroDivisionError:
        #    IPython.embed()
        ################################ WILD TESTING CODE HERE

        # Overkill sanity checking. This pipeline is so blazing fast, let's be paranoid.
        for param in ['start', 'stop', 'chromosome', 'strand']:
            assert(updated_dict[control_id][param] == control_data[param])

    return updated_dict


def main():
    parser = argparse.ArgumentParser(description='RIESLING: Rapid Identification of Enhancer Sites LInked to Nearby Genes.',
                                     epilog="Written by Nick Semenkovich <semenko@alum.mit.edu> for the Gordon Lab at "
                                            "Washington University in St. Louis: https://gordonlab.wustl.edu.",
                                     usage='%(prog)s [options]')

    parser.add_argument('--bam', dest="input_bam", metavar='input.bam', type=str,
                        help='The raw, original bam (e.g. the input to MACS).', required=True)

    parser.add_argument('--bed', dest="peak_bed", metavar='peaks.bed', type=str,
                        help='An input bed of peaks (e.g. MACS output).', required=True)

    parser.add_argument('--mask', dest="mask_bed", metavar='mask.bed', type=str,
                        help='Sites to mask from the .bed peaks, also .bed formatted.', required=False)

    parser.add_argument('--tss-window', dest="tss_window", metavar='tss_window', type=int,
                        help='Window (in bp) to mask around TSSes.', required=True)

    parser.add_argument('--stitch-window', dest="stitch_window", metavar='stitch_window', type=int,
                        help='Window (in bp) to stitch adjacent peaks together.', required=True)

    parser.add_argument('--force-dynstitch-window', dest="force_dynstitch_window", metavar='force_dynstitch_window', type=int,
                        help='Override the dynamic stitching algorithm with a fixed window (in bp)', required=False)

    parser.add_argument('--control-bam', dest="control_bam", metavar='control.bam', type=str,
                        help='A control bam with background/comparison data.', required=False)

    parser.add_argument('--out', '-o', dest="output", metavar='output_dir', type=str,
                        help='Output directory.', required=True)


    genome_choices = sorted(CONFIG['bowtie2_genomes'].keys())
    parser.add_argument('--genome', '-g', dest="genome", metavar='genome', type=str,
                        choices=genome_choices, help='RefSeq/TSS genome, one of: %s' % (', '.join(genome_choices)), required=True)

    other_ops = parser.add_argument_group('other options')

    other_ops.add_argument("--no-annotation", "-na", dest="no_annotation", default=False, action='store_true',
                           help="Skip final .bed annotation of nearby genes.")

    other_ops.add_argument("--no-bed-output", "-nb", dest="no_bed_output", default=False, action='store_true',
                           help="Skip .bed output creation. Implies --no-annotation")

    other_ops.add_argument("--verbose", "-v", dest="verbose", default=False, action='store_true')

    other_ops.add_argument("--no-log", "-nl", dest="nolog", default=False, action='store_true',
                           help="Do not create a log file.")

    other_ops.add_argument("--development-cache", dest="development_cache", default=False, action='store_true',
                           help="Cache intermediate files aggressively. Do *not* use in production.")


    args = parser.parse_args()

    # Convert input (either .bam or a path) to a list, with sanity checks (can we read? is list >0?)
    input_file_list = _script_helpers.validate_input_files(args.input_bam, mask=".bam")
    assert(len(input_file_list) == 1)

    output_path = _script_helpers.setup_output_path(args.output)

    ## Establish logging.
    log_main = _logshim.startLogger(verbose=args.verbose, noFileLog=args.nolog, outPath=output_path)

    # TODO: We've only been doing one-off runs -- can we even handle more than one input file?
    log_main.info("** Executing pipeline on %s input BAM file(s)." % (len(input_file_list)))
    log_main.info("** Writing to: %s/" % (output_path))

    # You can override the temp directory in .config.yaml.
    # Currently, this is a NOP, as the tempdir for pybedtools *is* the os tempdir.
    os_tmp_dir = tempfile.gettempdir()
    pybedtools.set_tempdir(os_tmp_dir)

    custom_tmp_dir = CONFIG['general']['custom_tmp_dir']
    if custom_tmp_dir is not None:
        if os.path.isdir(custom_tmp_dir) and os.access(custom_tmp_dir, os.W_OK):
            log_main.info('Using custom tmp directory: %s' % (custom_tmp_dir))
            pybedtools.set_tempdir(custom_tmp_dir)
        else:
            log_main.warn('Custom temp directory inaccessible, defaulting to OS.')


    # If we're tweaking things in development, use this cache structure.
    # This does *not* respect changes in underlying .bed/.bam files, or some calling flags (!!).
    # Be *very careful* using this flag -- it's most useful for tweaking graphs and TE/SE prediction code.

    bamliquidator_results = {}
    analysis_types = []

    if args.development_cache:
        log_main.warn('Trying to load results from development cache. Use *only* for development!')
        try:
            with open(output_path + '/global_pickle.cache', 'rb') as cache_fh:
                bamliquidator_results, analysis_types = cPickle.load(cache_fh)
                log_main.warn('Cache loaded successfully. DO NOT use this for final results!')
        except (IOError, cPickle.UnpicklingError):
            log_main.warn('  Cache corrupt or doesn\'t yet exist. Ignoring.')


    # If we don't have any results yet, generate them.
    if bamliquidator_results == {} or analysis_types == []:
        # Start by loading the UCSC refseq data for our genome of interest
        ucsc_refseq_bedtool, ucsc_name_mapping = load_refseq_ucsc_table(args.genome)

        # Handle peak masking, if we were asked to do so.
        # TODO: Integrate hot-spot identification here, too?
        # The goal here is to remove extreme density hotspots from input data.
        #  - You could use the ENCODE blacklists if available for your organism / assay
        #  - You could generate your own:
        #    * Hotspot masks can be done with a *folder* of BAMs, using bamliquidator_batch.py and a binsize of ~10k.
        #      Then run the flattener, e.g. https://github.com/BradnerLab/pipeline/wiki/bamliquidator#flattener
        #    * Rank order *density* from bamliquidator and exclude an insane IQR limit (e.g. >25)
        peak_bedtools = {}
        peak_bedtools['unstitched_with_tss'] = mask_bed_file(args.peak_bed, args.mask_bed)

        # Remove peaks within the TSS window
        peak_bedtools['unstitched'] = remove_transcription_start_sites(peak_bedtools['unstitched_with_tss'], ucsc_refseq_bedtool, ucsc_name_mapping, args.tss_window)

        # Merge adjacent features, based on a stitching window
        # TODO: Consider re-adding distinct gene boundaries that are overlapped?
        log_main.info('Performing both dynamic and static stitching.')

        dynstitched_unfiltered = stitch_bed_using_window(peak_bedtools['unstitched'], output_path, stitch_override=args.force_dynstitch_window)

        fixstitched_unfiltered = stitch_bed_using_window(peak_bedtools['unstitched'], output_path, stitch_override=args.stitch_window)


        # Now remove any TSSes re-introduced from the dynstitched and fixstitched populations.
        peak_bedtools['dynstitched'] = remove_transcription_start_sites(dynstitched_unfiltered, ucsc_refseq_bedtool, ucsc_name_mapping, args.tss_window, remove_entire_feature=False)

        peak_bedtools['fixstitched'] = remove_transcription_start_sites(fixstitched_unfiltered, ucsc_refseq_bedtool, ucsc_name_mapping, args.tss_window, remove_entire_feature=False)


        # Collect signal intensities using bamliquidator (this computes read density)
        log_main.info('Computing BAM read depth using un-/dyn-/fix-stitched peaks')
        analysis_types = ['unstitched', 'dynstitched', 'fixstitched']

        for type in analysis_types:
            bamliquidator_results[type] = compute_read_density(peak_bedtools[type], args.input_bam,
                                                               output_path, path_suffix='-%s' % type, development_cache=args.development_cache)

        if args.control_bam:
            log_main.info('** Control BAM present. Computing read density of groups on control BAM.')
            for type in analysis_types:
                bamliquidator_results['control_%s' % type] = compute_read_density(peak_bedtools[type], args.control_bam,
                                                                                  output_path, path_suffix='-%s-control' % type, development_cache=args.development_cache)

        if args.control_bam:
            log_main.info("** Removing background signal using control BED.")
            log_main.info("  You will have negative numbers in .bed output (!)")
            for type in analysis_types:
                bamliquidator_results['adjusted_%s' % type] = update_signal_using_control(bamliquidator_results[type], bamliquidator_results['control_%s' % type])

            # Expand the analyses to include these adjusted .bam signals, too.
            analysis_types += ['adjusted_%s' % elt for elt in analysis_types]

        # If requested, save the entire results object to our cache.
        if args.development_cache:
            log_main.warn('Writing to development cache! *Only* use during development!')
            with open(output_path + '/global_pickle.cache', 'wb') as cache_fh:
                cPickle.dump((bamliquidator_results, analysis_types), cache_fh, cPickle.HIGHEST_PROTOCOL)


    ## Now, we need to identify & segregate super enhancers from traditional enhancers
    log_main.info('Identifying super-enhancers (SEs) vs traditional enhancers (TEs)')

    SE_peaks = {}
    TE_peaks = {}
    # For control BAM adjusted negative values
    neg_SE_peaks = {}
    neg_TE_peaks = {}
    for type in analysis_types:
        log_main.info('Predicting enhancers for %s' % (type))

        full_output_path = '%s/%s' % (output_path, type)
        if not os.path.isdir(full_output_path):
            os.mkdir(full_output_path)

        if type.startswith('adjusted_'):
            SE_peaks[type], TE_peaks[type], neg_SE_peaks[type], neg_TE_peaks[type] = call_enhancers_with_control(bamliquidator_results[type], full_output_path, '%s' % type)
        else:
            SE_peaks[type], TE_peaks[type] = call_enhancers(bamliquidator_results[type], full_output_path, '%s' % type)


    if not args.no_bed_output:
        log_main.info('Exporting all data as .bed files')
        # Time to convert the SEs and TEs to output .bed files
        all_bed_outputs = {}
        pos_all_bed_outputs = {}
        neg_all_bed_outputs = {}

        se_bed_outputs = {}
        te_bed_outputs = {}

        neg_se_bed_outputs = {}
        neg_te_bed_outputs = {}

        for type in analysis_types:
            all_bed_outputs[type] = peak_id_to_bed(bamliquidator_results[type]).sort()

            se_bed_outputs[type] = peak_id_to_bed(bamliquidator_results[type], peak_ids=SE_peaks[type]).sort()
            te_bed_outputs[type] = peak_id_to_bed(bamliquidator_results[type], peak_ids=TE_peaks[type]).sort()

            if type.startswith('adjusted_'):
                pos_all_bed_outputs[type] = peak_id_to_bed(bamliquidator_results[type], peak_ids=SE_peaks[type] + TE_peaks[type]).sort()
                neg_all_bed_outputs[type] = peak_id_to_bed(bamliquidator_results[type], peak_ids=neg_SE_peaks[type] + neg_TE_peaks[type]).sort()

                neg_se_bed_outputs[type] = peak_id_to_bed(bamliquidator_results[type], peak_ids=neg_SE_peaks[type]).sort()
                neg_te_bed_outputs[type] = peak_id_to_bed(bamliquidator_results[type], peak_ids=neg_TE_peaks[type]).sort()

        # Save those bed outputs
        for type in analysis_types:
            full_output_path = '%s/%s' % (output_path, type)
            all_bed_outputs[type].saveas('%s/%s_all.bed' % (full_output_path, type))

            se_bed_outputs[type].saveas('%s/%s_SEs_pos.bed' % (full_output_path, type))
            te_bed_outputs[type].saveas('%s/%s_TEs_pos.bed' % (full_output_path, type))

            if type.startswith('adjusted_'):
                pos_all_bed_outputs[type].saveas('%s/%s_all_pos.bed' % (full_output_path, type))
                neg_all_bed_outputs[type].saveas('%s/%s_all_neg.bed' % (full_output_path, type))

                neg_se_bed_outputs[type].saveas('%s/%s_SEs_neg.bed' % (full_output_path, type))
                neg_te_bed_outputs[type].saveas('%s/%s_TEs_neg.bed' % (full_output_path, type))


        if not args.no_annotation:
            # And now let's sort & annotate all these files, too.
            genome_gff = CONFIG['gffs'][args.genome]

            log_main.info('Finding nearby genes for .annotated.bed output (this takes a minute)')
            for type in analysis_types:
                full_output_path = '%s/%s' % (output_path, type)
                log_main.info('  Working on: %s' % type)
                # D=ref reports distance, with 'upstream' negative
                # t=first reports the first entry from the reference .gtf if multiple match (otherwise report all)
                all_bed_outputs[type].closest(genome_gff, D='ref', t='first').saveas('%s/%s_all.annotated.bed' % (full_output_path, type))

                # Repeating unnecessary work. TODO: Prune/disable/reduce this?
                se_bed_outputs[type].closest(genome_gff, D='ref', t='first').saveas('%s/%s_SEs_pos.annotated.bed' % (full_output_path, type))
                te_bed_outputs[type].closest(genome_gff, D='ref', t='first').saveas('%s/%s_TEs_pos.annotated.bed' % (full_output_path, type))

                if type.startswith('adjusted_'):
                    # TODO: Skip all this repeat stuff?
                    pos_all_bed_outputs[type].closest(genome_gff, D='ref', t='first').saveas('%s/%s_all_pos.annotated.bed' % (full_output_path, type))
                    neg_all_bed_outputs[type].closest(genome_gff, D='ref', t='first').saveas('%s/%s_all_neg.annotated.bed' % (full_output_path, type))

                    neg_se_bed_outputs[type].closest(genome_gff, D='ref', t='first').saveas('%s/%s_SEs_neg.annotated.bed' % (full_output_path, type))
                    neg_te_bed_outputs[type].closest(genome_gff, D='ref', t='first').saveas('%s/%s_TEs_neg.annotated.bed' % (full_output_path, type))


        else:
            log_main.warning('Skipping annotation as requested.')

    else:
        log_main.warning('Skipping bed creation as requested.')


    # Cleanup this session's bedtool objects.
    pybedtools.cleanup(verbose=False, remove_all=False)


if __name__ == '__main__':
    main()
