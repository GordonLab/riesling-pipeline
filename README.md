## The RIESLING ATAC-seq Pipeline [![Build Status](https://travis-ci.org/GordonLab/riesling-pipeline.svg?branch=master)](https://travis-ci.org/GordonLab/riesling-pipeline)

The RIESLING (Rapid Identification of EnhancerS LInked to Nearby Genes) ATAC-seq pipeline is designed to be a simplistic set of standalone scripts for analyzing ATAC-seq data.

You may find it particularly useful for identifying and stratifying super-enhancers, though can also be leveraged for differential accessibility analysis using DESeq2.

Since this was originally developed in 2014-2015, a number of other packages have been developed, notably the [Kundaje lab pipeline](https://github.com/kundajelab/atac_dnase_pipelines) and the [Ren lab single-cell ATAC pipeline](https://github.com/r3fang/scATAC). The Kundaje lab pipeline in particular includes other features (e.g. IDR analysis) which you may find more useful if you aren't interested in super-enhancers / enhancer clusters.

================

## Getting started
1. Clone this repo: `git clone https://github.com/GordonLab/riesling-pipeline.git`
2. `cd riesling-pipeline`
3. Install the Python dependencies: `pip install --user -U -r requirements.txt`
