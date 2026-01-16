# =============================================================================
# CysRichPepScanR - NCR-like Peptide Discovery Pipeline (Main Script)
# 
# This script performs motif-based and homology-based screening of cysteine-rich
# peptides from predicted proteins (GeMoMa or ORFipy output).
# Optimized for NCR peptides in IRLC legumes; adaptable to other CRP families.
#
# Author: Sneha Bidvi (20micro.sneha@iscm.ac.in)
# Based on Bidvi et al. (2026) - A bioinformatics pipeline for Screening Nodule-specific Cysteine-Rich (NCR) like peptides from Trigonella foenum-graecum and Medicago truncatula genomes
# Date: January 2026
# =============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# 1. USER-DEFINED PARAMETERS (customize here)
# ──────────────────────────────────────────────────────────────────────────────


setwd("E:\\trial")


input_proteins_fasta   <- "predicted_proteins.fasta"           # From GeMoMa/ORFipy
output_dir             <- "results"                            # Create if missing
ncr_reference_faa      <- "databases/NCR_peptides.faa"         # Known NCR sequences
blast_db_name          <- "databases/NCRdb"                   # Will be created
threads                <- 16                                   # Adjust to your CPU

min_length             <- 20
max_length             <- 180
min_cysteines          <- 4

motifs <- c(
  "CX{5}CX{0,41}CX{0,16}CX{4}CX{1}C",   # 6-cysteine pattern
  "CX{5}CX{1,33}CX{4}C"                 # 4-cysteine pattern
)

prefix_6c <- "ncr6c_"
prefix_4c <- "ncr4c_"

# SignalP options
signalp_organism <- "eukarya"
signalp_mode     <- "fast"   # or "slow" for more accurate borders

# ──────────────────────────────────────────────────────────────────────────────
# 2. SETUP & PACKAGES
# ──────────────────────────────────────────────────────────────────────────────

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(Biostrings)
  library(ampir)      # optional, but used in original
  library(dplyr)
  library(readr)      # for safer CSV reading
})

# Helper: check file exists
check_file <- function(f, desc = basename(f)) {
  if (!file.exists(f)) stop(sprintf("File not found: %s (%s)", f, desc))
}

check_file(input_proteins_fasta, "input proteins FASTA")
check_file(ncr_reference_faa,    "NCR reference database")

setwd(output_dir)  # All outputs go here

# ──────────────────────────────────────────────────────────────────────────────
# 3. FILTER PEPTIDES: length 20-180 aa + ≥4 cysteines
# ──────────────────────────────────────────────────────────────────────────────

message("Reading predicted proteins...")
peptides <- readAAStringSet(input_proteins_fasta)

pep_lengths <- width(peptides)
cys_counts  <- letterFrequency(peptides, "C")

message(sprintf("Total peptides before filter: %d", length(peptides)))

peptides_filt <- peptides[
  pep_lengths >= min_length &
    pep_lengths <= max_length &
    cys_counts  >= min_cysteines
]

message(sprintf("Peptides after filter: %d", length(peptides_filt)))

filtered_fasta <- "filtered_peptides_20-180aa_>=4Cys.fasta"
writeXStringSet(peptides_filt, filtered_fasta)

# Use unique sequences for motif screening
seqs_unique <- unique(as.character(peptides_filt))

# ──────────────────────────────────────────────────────────────────────────────
# 4. MOTIF-BASED SCREENING (4C & 6C patterns)
# ──────────────────────────────────────────────────────────────────────────────

screen_peptides <- function(seqs, motifs) {
  matches <- list()
  for (mot in motifs) {
    pattern <- gsub("X", ".", mot)          # regex-friendly
    hits    <- grep(pattern, seqs, value = TRUE)
    matches[[mot]] <- hits
    message(sprintf("Motif '%s' matched %d sequences", mot, length(hits)))
  }
  matches
}

message("Screening motifs...")
motif_matches <- screen_peptides(seqs_unique, motifs)

# Write motif-specific FASTA files
for (i in seq_along(motifs)) {
  mot <- motifs[i]
  hits <- motif_matches[[mot]]
  if (length(hits) == 0) next
  
  prefix <- if (mot == motifs[1]) prefix_6c else prefix_4c
  fname  <- paste0("motif_hits_", gsub("[^A-Za-z0-9]", "_", mot), ".fasta")
  
  con <- file(fname, "w")
  for (j in seq_along(hits)) {
    cat(sprintf(">%s%d\n%s\n", prefix, j, hits[j]), file = con)
  }
  close(con)
  message(sprintf("Wrote %d sequences to %s", length(hits), fname))
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. DEDUPLICATE 4C against 6C (avoid overlap)
# ──────────────────────────────────────────────────────────────────────────────

fasta_6c <- paste0("motif_hits_", gsub("[^A-Za-z0-9]", "_", motifs[1]), ".fasta")
fasta_4c <- paste0("motif_hits_", gsub("[^A-Za-z0-9]", "_", motifs[2]), ".fasta")

check_file(fasta_6c); check_file(fasta_4c)

seqs6c <- readAAStringSet(fasta_6c)
seqs4c <- readAAStringSet(fasta_4c)

seqstr_6c <- as.character(seqs6c)
seqstr_4c <- as.character(seqs4c)

seqs4c_dedup <- seqs4c[!seqstr_4c %in% seqstr_6c]

writeXStringSet(seqs4c_dedup, "motif_hits_4c_dedup.fasta")
message(sprintf("4C after deduplication vs 6C: %d", length(seqs4c_dedup)))

# ──────────────────────────────────────────────────────────────────────────────
# 6. HOMOLOGY SEARCH with BLASTp
# ──────────────────────────────────────────────────────────────────────────────

run_cmd <- function(cmd, args, desc) {
  message(sprintf("Running: %s %s", cmd, paste(args, collapse=" ")))
  status <- system2(cmd, args)
  if (status != 0) stop(sprintf("Command failed: %s (%s)", desc, status))
}

# Build BLAST DB (only if missing)
if (!file.exists(paste0(blast_db_name, ".pdb"))) {
  message("Building BLAST database...")
  run_cmd("makeblastdb", c(
    "-in",  shQuote(ncr_reference_faa),
    "-dbtype", "prot",
    "-out", blast_db_name
  ), "makeblastdb")
}

# Run BLASTp for both sets
for (cys_type in c("6c", "4c")) {
  query_fasta <- if (cys_type == "6c") fasta_6c else "motif_hits_4c_dedup.fasta"
  out_tsv     <- sprintf("blastp_hits_%s.tsv", cys_type)
  
  check_file(query_fasta)
  
  run_cmd("blastp", c(
    "-query",   shQuote(query_fasta),
    "-db",      blast_db_name,
    "-out",     out_tsv,
    "-evalue",  "1e-3",
    "-outfmt",  "6",
    "-num_threads", as.character(threads)
  ), sprintf("BLASTp %s", cys_type))
  
  # Convert to CSV with headers
  dt <- read.delim(out_tsv, header = FALSE, stringsAsFactors = FALSE)
  colnames(dt) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore")
  write_csv(dt, sprintf("blastp_hits_%s.csv", cys_type))
}

# ──────────────────────────────────────────────────────────────────────────────
# 7. FILTER BLAST HITS (bitscore >=50, best per query by evalue)
# ──────────────────────────────────────────────────────────────────────────────

filter_blast <- function(csv_file) {
  dt <- read_csv(csv_file, show_col_types = FALSE)
  dt <- dt %>% filter(bitscore >= 50)
  top <- dt %>%
    group_by(qseqid) %>%
    slice_min(evalue, with_ties = FALSE) %>%
    ungroup()
  top$qseqid
}

ids_6c <- filter_blast("blastp_hits_6c.csv")
ids_4c <- filter_blast("blastp_hits_4c.csv")

writeLines(ids_6c, "blastp_top_ids_6c.txt")
writeLines(ids_4c, "blastp_top_ids_4c.txt")

# Extract candidate sequences
run_cmd("seqkit", c("grep", "-f", "blastp_top_ids_6c.txt", shQuote(fasta_6c),
                    ">", "candidates_6c.fasta"), "seqkit 6c")
run_cmd("seqkit", c("grep", "-f", "blastp_top_ids_4c.txt", "motif_hits_4c_dedup.fasta",
                    ">", "candidates_4c.fasta"), "seqkit 4c")

# ──────────────────────────────────────────────────────────────────────────────
# 8. SIGNALP 6 PREDICTION (mature peptides)
# ──────────────────────────────────────────────────────────────────────────────

for (cys_type in c("6c", "4c")) {
  infile  <- sprintf("candidates_%s.fasta", cys_type)
  outdir  <- sprintf("signalp_%s", cys_type)
  
  check_file(infile)
  
  run_cmd("signalp6", c(
    "--fastafile", shQuote(infile),
    "--organism", signalp_organism,
    "--output_dir", outdir,
    "--format", "txt",
    "--mode", signalp_mode
  ), sprintf("SignalP %s", cys_type))
}

# ──────────────────────────────────────────────────────────────────────────────
# 9. COMBINE MATURE PEPTIDES & FINAL OUTPUT
# ──────────────────────────────────────────────────────────────────────────────

mature_6c <- readAAStringSet(file.path("signalp_6c", "processed_entries.fasta"))
mature_4c <- readAAStringSet(file.path("signalp_4c", "processed_entries.fasta"))

all_mature <- c(mature_6c, mature_4c)
all_mature <- all_mature[!duplicated(as.character(all_mature))]

# Number them
names(all_mature) <- paste0("ncr_", seq_along(all_mature))

final_fasta <- "final_ncr_like_peptides.fasta"
writeXStringSet(all_mature, final_fasta)

message(sprintf("Final NCR-like candidates: %d sequences", length(all_mature)))
message(sprintf("Saved to: %s", final_fasta))

# ──────────────────────────────────────────────────────────────────────────────
# End of pipeline
# ──────────────────────────────────────────────────────────────────────────────