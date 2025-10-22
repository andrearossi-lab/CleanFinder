CleanFinder

CleanFinder: A Browser-Native Suite for Genome Editing Analysis

CleanFinder is a lightweight, all-in-one web application for designing, validating, and analyzing genome editing experiments. It runs 100% in your web browser, requiring no installation, no server, and no data upload. This client-side approach guarantees that your sensitive genomic data never leaves your computer.

CleanFinder is a comprehensive suite that integrates multiple tools into a single, seamless workflow.

Core Modules

1. gDNA Finder (Amplicon Finder)
This module retrieves genomic data in real-time from online repositories (like Ensembl), ensuring you always use the most up-to-date reference genome.
Paired-Primer & Single/sgRNA Modes: Retrieve amplicons using either primer pairs or a single guide sequence.
Genomic Context: Automatically identifies and lists genes within the retrieved amplicon region.
Allele Dropout Prevention: By visualizing your primers relative to the cut site, this module helps you design robust experiments that prevent allele dropout, a common artifact where large deletions remove a primer binding site.

2. Analyze FASTQ
The core analysis engine for quantifying your editing results from short-read sequencing (e.g., Illumina).
Robust Anchor-Based Analysis: Uses a buffer-based anchoring system to reliably identify on-target reads.
Functional Classification: Automatically classifies every read into intuitive, functionally-relevant categories:
    - Wild-type, HDR, Substitution/SNP
    - In-frame and Out-of-frame Indels
Rich Visualizations: Generates interactive graphs and tables for each sample, including:
    - Allele frequency charts (grouped by functional type or by specific alleles)
    - Indel size distribution plots
    - A Positional Mutation Profile heatmap.
HTML Report: Download a complete, self-contained HTML report of your analysis for sharing and archiving.

3. Allelic Dropout SNP Analyzer
A dedicated module to detect allele dropout in your experiments. This tool uses long-read FASTQ data (e.g., PacBio/Nanopore) from a heterozygous sample to find a high-confidence heterozygous SNP. This SNP can then be used as a marker to confirm if one allele is "dropping out" (i.e., failing to amplify) in your separate short-read editing analysis.

4. Transcript Finder
An in-silico PCR tool for cDNA. It checks your primers against all known transcripts (isoforms) of a gene to predict amplicon sizes. This is ideal for validating qPCR primer specificity and ensuring you are amplifying the correct isoforms.

5. Genome PCR
An in-silico PCR tool that blasts your primer pair against the entire human genome to check for specificity and predict potential off-target amplicons.

6. Rev/RevComp Tool
A simple utility to quickly calculate the reverse and reverse-complement of any DNA sequence.

---

Alignment Engine

CleanFinder's analysis modules (FASTQ Analyzer and Allelic Dropout Analyzer) are powered by a client-side JavaScript implementation of the Needleman–Wunsch global alignment algorithm.

Algorithm: Needleman-Wunsch with affine gap penalties.
- Dynamic Scoring Matrix: The scoring is dynamically adjusted based on the selected enzyme for more biologically accurate results:
  - SpCas9: { match: +5, mismatch: -4, gap open: -10, gap extend: -1 }
  - Cas12a: { match: +4, mismatch: -3, gap open: -8, gap extend: -1 }
  - Default: { match: +2, mismatch: -1, gap open: -5, gap extend: -1 }

This algorithm performs a highly accurate, per-read alignment between the reference sequence and each sequencing read, allowing it to precisely characterize indels near the CRISPR cut site.

Availability

CleanFinder is distributed as a single HTML file with embedded JavaScript and can be run in any modern web browser without installation. The live application is accessible by visiting andrearossi.de and launching the tool via the prompt open cleanfinder. The core JavaScript alignment algorithm is available under the MIT license at https://github.com/andrearossi-lab.
