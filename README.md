## The RIESLING ATAC-seq Pipeline [![Build Status](https://travis-ci.org/GordonLab/riesling-pipeline.svg?branch=master)](https://travis-ci.org/GordonLab/riesling-pipeline)

The RIESLING (Rapid Identification of EnhancerS LInked to Nearby Genes) ATAC-seq pipeline is designed to be an efficient set of standalone scripts for quickly analyzing ATAC-seq data.

You may find it particularly useful for identifying and stratifying super-enhancers, though can also be leveraged for differential accessibility analysis using DESeq2.

Since this was originally developed in 2014-2015, a number of other packages have been developed, notably the [Kundaje lab pipeline](https://github.com/kundajelab/atac_dnase_pipelines) and the [Ren lab single-cell ATAC pipeline](https://github.com/r3fang/scATAC). The Kundaje lab pipeline in particular includes other features (e.g. IDR analysis) which you may find more useful if you aren't interested in super-enhancers / enhancer clusters.

================

## Getting started
1. Clone this repo: `git clone https://github.com/GordonLab/riesling-pipeline.git`
2. `cd riesling-pipeline`
3. Install the Python dependencies: `pip install --user -U -r requirements.txt`

=================
/*
# A Working Example

** (Coming soon) **
*/

## Expected Inputs & Pre-processing Data

The heart of RIESLING and related code expects:

* A .bam file of your aligned ATAC-seq data (e.g. from bowtie2)
* A .bed of peaks (from MACS, HOMER, or any other standard peak caller)

There are lots of ways to pre-process data -- please use whatever approach you prefer. A basic set of standalone
scripts to help pre-process paired-end sequencing data is provided, detailed below under 'Pre-processing Tools'.

/*
## Calling Super-Enhancers, Stretch Enhancers, and more

** (Full example coming soon) **

## Differential Accessibility Analyses (e.g. with DESeq2)

** (Full example coming soon) **
*/

## Pre-processing Tools

There are lots of valid way to pre-process data. These standalone scripts may be of use, but please pay careful attention
to the peak calling settings you use -- as the defaults may not be applicable to your experimental approach.


* (Optional) ./0-merge-fastq.py: Intelligently merge across lanes, for multiple-lane samples (e.g. one multiplexed sample loaded into multiple lanes).
This will only be useful if you've loaded the *same, multiplexed sample* into multiple lanes of an Illumina flowcell.

This script concatenates the same index's pared-end files (L*_R*_* .fastq.gz) across multiple lanes into one set of PE files per-sample.


* ./1-map-to-genome.py: Run bowtie2 across a folder of paired-end sequence data.

* ./2-sanitize-bam.py: Clean up .bam files, including ATAC-specific fixes.
This includes: chrM removal, quality filtering, ENCODE blacklist removal, and PCR duplicate removal

* (Optional) ./3-OPTIONAL-merge-bam-rmdup.py: A helper script to blindly concatenate and deduplicate multiple sets of BAMs.

* ./4-call-peaks.py: Run both macs14 and macs2 on .bam files.

By default, this runs both macs14 and macs2, and operates on directories of .bam files.
