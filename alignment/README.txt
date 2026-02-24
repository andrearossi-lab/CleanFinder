Synthetic Illumina SE datasets (5,000 reads each), 10% edited, based on reference_amplicon_200bp.fasta

>amplicon_ref_200bp_with_gRNA_PAM
AAGAAGGTCAATGTCCAATCTAGAACTCCAAACGAGTGTCCGCTTGAAGTTCAATTCGTAATAGATCTGACACACATTCG
GAAGGATAGCTGACCTGATCGTACGTTGAGGTCACCGACAGACGGGACCACCCCGAACGGAAGATTATCCGGGGATCTAT
AAGAGATCACAGTCTAGTAGGAACAAAACTAGGACGGCTC

gRNA (20nt): GCTGACCTGATCGTACGTTG
PAM: AGG
Target start (0-based): 88
Cas9 cut position (conceptual): between 104 and 105

Datasets:
 - sim_SE200_5k_KI15_10pct.fastq: 15 bp insertion at cut site (KI15). 
Inserted sequence: GCTAGTCCGATGACG

Sequence for KI
GTACGgctagtccgatgacgTTGAG


 - sim_SE200_5k_CBE_CtoT_10pct.fastq: C->T at reference index 93 (protospacer window 4-8)

GCTGAtCTGAT. Window open 6 closes at 6.

 - sim_SE200_5k_ABE_AtoG_10pct.fastq: A->G at reference index 92 (protospacer window 4-8)

AGCTGgCCTGA. Window open 6 closes at 6.



Notes:
 - WT reads are 200 bp; KI reads are 215 bp.
 - Quality strings are constant 'I' (Phred 40)


PE: Edits (0-based index, ref>alt): 106:T>C, 110:G>T, 111:T>A


TACGTCGAGTACACCG
