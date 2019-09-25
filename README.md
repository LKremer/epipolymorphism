# calc_epipolymorphism.py
A command-line tool to estimate epigenetic polymorphism from bisulfite-converted sequencing reads.

Epipolymorphism is calculated according to the formula in Landan et. al 2012 (doi: 10.1038/ng.2442).

## Requirements
requires Python 3 and the Python package `pysam`.

## Output
The output is a tab-separeted text file (tsv) containing the following columns:
`chr, start, end, epipolymorphism, mfrac, coverage, spanning_reads`
Where `chr`, `start` and `end` are the coordinates of the CpG-group (by default 4 CpGs),
`epipolymorphism` is the estimated epipolymorphism, `mfrac` is the fraction of methylated
reads in the CpG-group (only using reads that span the whole region as when computing
epipolymorphism), `coverage` is the total coverage of the region including non-spanning reads,
and `spanning_reads` is the number of reads spanning the region (used for calculating `mfrac`
and `epipolymorphism`).

## Usage

Example:

    python3 calc_epipolymorphism.py bs_reads.bam cpg_locations.bed epipoly_output.tsv

Use `./calc_epipolymorphism.py --help` to see the usage instructions:

    usage: calc_epipolymorphism.py [-h] [--check-fasta CHECK_FASTA]
                                   [--min-cluster-width MIN_CLUSTER_WIDTH]
                                   [--max-cluster-width MAX_CLUSTER_WIDTH]
                                   [--cluster-n CLUSTER_N]
                                   bam_file bed_file output_tsv

    Calculates epipolymorphism as defined in Landan et al. 2012 from a SAM or BAM
    file of bisulfite-converted DNA reads

    positional arguments:
      bam_file              BAM or SAM file of bisulfite-converted, aligned reads
      bed_file              bed file of genomic CpG sites of interest (usually all
                            CpG sites)
      output_tsv            output tsv file path

    optional arguments:
      -h, --help            show this help message and exit
      --check-fasta CHECK_FASTA
                            optionally, use the genome fasta to assert that the
                            CpG locations in the bed file are really CpG sites
      --min-cluster-width MIN_CLUSTER_WIDTH
                            CpG clusters below this width will not be reported
      --max-cluster-width MAX_CLUSTER_WIDTH
                            CpG clusters above this width will not be reported
                            (default: 100; this should be less or equal to the
                            read length)
      --cluster-n CLUSTER_N
                            number of CpG sites in a cluster (default: 4; this is
                            the usual value)
