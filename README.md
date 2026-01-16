# CysRichPepScanR

A modular, motif- and homology-based pipeline for discovery of cysteine-rich secreted peptides (CRPs) from genome or transcriptome data. The workflow is optimized for NCR peptides but is easily extensible to defensins, thionins, cyclotides, and other CRP families via parameter and motif customization.

This repository contains an R-based bioinformatics workflow for identifying Nodule-specific Cysteine-Rich (NCR)-like peptides from legume genomes or transcriptomes. It uses motif-based and homology-based screening, as described in the paper "R-Guided Discovery of NCR-Like Peptides in Trigonella foenum-graecum and Medicago truncatula Genomes Using Motif and Homology-Based Screening" (Bidvi et al., 2026).

The pipeline predicts open reading frames (ORFs) or genes, filters peptides by length (20-180 aa) and cysteine content (>=4), screens for conserved cysteine motifs (4C or 6C patterns), performs BLASTp against a custom NCR database, filters hits, predicts signal peptides with SignalP, and outputs final NCR-like sequences.
While optimized for NCR peptides in IRLC legumes, components like ORF prediction, length/cysteine filtering, motif screening, and signal peptide prediction can be adapted for other cysteine-rich peptide families (e.g., defensins, thionins, cyclotides) by tuning parameters, motifs, or reference sequences.

**Prerequisites**

Operating System: Linux/Unix (tested on Ubuntu; paths assume absolute Linux paths).
R: Version 4.0+ with packages:
Biostrings (from Bioconductor: BiocManager::install("Biostrings")).
ampir (CRAN: install.packages("ampir")).
dplyr (CRAN: install.packages("dplyr")).
readxl (CRAN: install.packages("readxl")).

Conda/Miniconda: For managing environments and tools.

**External Tools:**

GeMoMa: For genome gene prediction (install via Conda: conda install -c bioconda gemoma).
ORFipy: For transcriptome ORF prediction (install via Conda: conda install -c bioconda orfipy).
BLAST+: For homology search (install via Conda: conda install -c bioconda blast).
Seqkit: For sequence extraction (install via Conda: conda install -c bioconda seqkit).
SignalP 6: For signal peptide prediction. This is a Python-based tool available for academic users (free with registration) or commercial licensing from DTU Health Tech. It predicts signal peptides and cleavage sites, essential for confirming secretory NCR-like peptides.

Custom NCR Database: A FASTA file of known NCR peptides (e.g., from Montiel et al., 2017; place in databases/NCR_peptides.faa).
Reference Genomes/Annotations: For GeMoMa (e.g., Medicago and Trifolium; paths defined in config.R).
Hardware: Multi-core CPU (uses 16 threads by default) and sufficient storage for large FASTA files.


**Citation**
If using this pipeline, please cite:

Bidvi S, Choure R, Jadhav R, Padul M, Mandavkar S, Bhadane A, and Posam M. A bioinformatics
pipeline for Screening Nodule-specific Cysteine-Rich (NCR) like peptides from Trigonella foenum-
graecum and Medicago truncatula genomes. In: (2025). doi: https://doi.org/10.21203/rs.3.rs-
7218712/v1

**Contact** 

For questions, contact Sneha Bidvi (20micro.sneha@iscm.ac.in)
