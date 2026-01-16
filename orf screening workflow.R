setwd("/home/rajendra/Documents")

# ===== LOAD CONFIG =====
source("config.R")  # Must define: data_type, input_fasta, output_dir, conda_bin_path, gemoma_outdir, ref_genomes, orfipy_path, gemoma_path

# ===== ENVIRONMENT SETUP =====
Sys.setenv(PATH = paste(conda_bin_path, Sys.getenv("PATH"), sep = ":"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===== ORFIPY FUNCTION =====
run_orfipy <- function(input_fasta, outdir, prefix = "predicted_orfs") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  args <- c(
    shQuote(input_fasta),
    "--outdir", outdir,
    "--pep", paste0(prefix, "_pep.fa"),
    "--dna", paste0(prefix, "_dna.fa"),
    "--min", "60",
    "--max", "5400"
  )
  system2(orfipy_path, args = args)
}

run_gemoma <- function() {
  dir.create(gemoma_outdir, recursive = TRUE, showWarnings = FALSE)  # Ensure output directory exists
  
  args <- c(
    "GeMoMaPipeline",
    "threads=16",
    paste0("outdir=", gemoma_outdir),
    "GeMoMa.Score=ReAlign",
    "AnnotationFinalizer.r=NO",
    "o=true",
    paste0("t=", input_fasta)
  )
  
  for (ref in ref_genomes) {
    args <- c(args,
              "s=own",
              paste0("i=", ref$id),
              paste0("a=", ref$annotation),
              paste0("g=", ref$genome))
  }
  
  system2(gemoma_path, args = args)
}

# ===== CONDITIONAL EXECUTION =====
if (data_type == "genome") {
  message("Genome input → Running GeMoMa...")
  run_gemoma()
} else if (data_type == "transcriptome") {
  message("Transcriptome input → Running ORFipy...")
  run_orfipy(input_fasta, file.path(output_dir, "orfipy_out"))
} else {
  stop("Invalid data_type. Must be 'genome' or 'transcriptome'")
}


# ===== PEPTIDE SCREENING FUNCTION =====
filter_peptides <- function(fasta_file, min_len = 20, max_len = 180, min_cys = 4, output_file = "filtered_peptides.fa") {
  library(Biostrings)
  
  aa_seqs <- readAAStringSet(fasta_file)
  
  # Filter by length and cysteine content
  valid_seqs <- aa_seqs[
    width(aa_seqs) >= min_len &
      width(aa_seqs) <= max_len &
      vcountPattern("C", aa_seqs) >= min_cys
  ]
  
  message(length(valid_seqs), " peptides passed the filter.")
  
  writeXStringSet(valid_seqs, filepath = output_file)
}

# ===== FILTER BASED ON SOURCE =====
if (data_type == "genome") {
  peptide_input_file <- file.path(gemoma_outdir, "predicted_proteins.fasta")
  output_filtered_file <- file.path(gemoma_outdir, "filtered_peptides.fa")
} else if (data_type == "transcriptome") {
  peptide_input_file <- file.path(output_dir, "orfipy_out", "predicted_orfs_pep.fa")
  output_filtered_file <- file.path(output_dir, "orfipy_out", "filtered_peptides.fa")
}

# Run filter
filter_peptides(peptide_input_file, output_file = output_filtered_file)

