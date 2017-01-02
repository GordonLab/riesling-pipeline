## Optional & Helper Scripts

* `0-merge-fastq.py`: Intelligently merge across lanes, for multiple-lane samples (e.g. one multiplexed sample loaded into multiple lanes).
This will only be useful if you've loaded the *same, multiplexed sample* into multiple lanes of an Illumina flowcell.

This script concatenates the same index's pared-end files (L*_R*_* .fastq.gz) across multiple lanes into one set of PE files per-sample.

* `3-merge-bam-rmdup.py`: A helper script to blindly concatenate and deduplicate multiple sets of BAMs.

* `hdf5_to_counts_table.py`: This script exists to convert .hdf5 files from bamliquidator into counts tables readable by R.
