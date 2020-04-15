# boll_weevil_pop_gen
Scripts used for QC, SNP calling and filtering, and analysis of boll weevil ddRADseq data.

This repository contains scripts for the analysis of ddRADseq data generated from boll weevil, Anthonomus grandis, for the purposes of population genomic analysis. Computationally intensive segments of these scripts were performed on the Texas A&M University High Performance Research Computing cluster "Ada." QC steps included here are filtering bacterial contaminants using Kraken, trimming using Trimmomatic, and removal of mitochondiral sequences using fastQscreen. SNP calling was performed with dDocent and SNP filtering was performed with vcftools. Analysis of the filtered dataset was carried out using R.
