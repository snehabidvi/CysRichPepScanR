# ===== BASIC SETUP =====
data_type   <- "genome"  # "genome" or "transcriptome"
input_fasta <- "/path/to/input_file.fasta"

output_dir    <- "/home/rajendra/Documents/results"
gemoma_outdir <- file.path(output_dir, "gemoma_out")

# ===== PATHS TO EXECUTABLES (absolute, adjust if needed) =====
gemoma_env_bin <- "/home/user/miniconda3/envs/gemoma_env/bin"
orfipy_path    <- "/home/user/orfipy-env/bin/orfipy"
gemoma_path    <- "/home/user/miniconda3/envs/gemoma_env/bin/GeMoMa"

# Try to detect conda base; fall back to a sensible default
detect_conda_base <- function() {
  b <- tryCatch(system2("bash", c("-lc", "conda info --base"), stdout = TRUE), error = function(e) "") 
  if (length(b) > 0 && nzchar(b[1])) b[1] else "/home/user/miniconda3"
}
conda_base     <- detect_conda_base()
conda_bin_path <- file.path(conda_base, "bin")  # <-- DEFINE BEFORE USING

# ===== SAFE HELPERS =====
prepend_path <- function(dir) {
  if (nzchar(dir) && dir.exists(dir)) {
    Sys.setenv(PATH = paste(dir, Sys.getenv("PATH"), sep = ":"))
  }
}
req_file <- function(f, what) {
  if (!file.exists(f)) stop(sprintf("Missing %s: %s", what, f))
}

# ===== PREPEND PATHS (order matters) =====
prepend_path(gemoma_env_bin)
prepend_path(conda_bin_path)

# ===== INPUT & OUTPUT CHECKS =====
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gemoma_outdir, showWarnings = FALSE, recursive = TRUE)
req_file(input_fasta, "input_fasta")

if (data_type == "transcriptome" && grepl("\\.(fa|fasta)(\\.gz)?$", input_fasta, ignore.case = TRUE) &&
    grepl("genome", basename(input_fasta), ignore.case = TRUE)) {
  warning("data_type='transcriptome' but input_fasta looks like a genome file. Double-check.")
}

# ===== REFERENCE GENOMES FOR GeMoMa =====
ref_genomes <- list(
  list(id = "ref1",
       annotation = "/path/to/ref1_annotation.gff3",
       genome     = "/path/to/ref1.fa"),
  list(id = "ref2",
       annotation = "/path/to/ref2_annotation.gff3",
       genome     = "/path/to/ref2.fa")
)

# ===== EXECUTABLE PRESENCE CHECKS (informative, not fatal) =====
exe <- Sys.which(c(GeMoMa = "GeMoMa", orfipy = "orfipy", mmseqs = "mmseqs"))
if (nzchar(orfipy_path) && !file.exists(orfipy_path) && !nzchar(exe["orfipy"])) {
  warning("orfipy not found at orfipy_path and not in PATH.")
}
if (nzchar(gemoma_path) && !file.exists(gemoma_path) && !nzchar(exe["GeMoMa"])) {
  warning("GeMoMa not found at gemoma_path and not in PATH.")
}
