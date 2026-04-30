#!/usr/bin/env Rscript

# Analyze one paired monomer/dimer input and write per-sequence topology fractions.
# The output `state` column refers specifically to assembly state: monomer or dimer.
# Package roles:
# - dplyr/readr: bind output rows and write CSV
# - network: construct the undirected contact network object expected by fibril_assay
# - orca: required indirectly because fibril_assay uses ORCA graphlet orbit counts

script_path <- normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- dirname(sub("^--file=", "", script_path))

# Minimal package set: data IO, row-binding, and undirected network construction.
required_packages <- c('dplyr', 'readr', 'network', 'orca')
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop('Missing packages: ', paste(missing, collapse = ', '))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(network)
})

source(file.path(script_dir, 'fibril_assay.R'))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop('Usage: Rscript topology_pair_analysis.R <monomer_rda> <dimer_single_node_rda> <label> <output_csv>')
}

monomer_file <- args[[1]]
dimer_file <- args[[2]]
label <- args[[3]]
out_csv <- args[[4]]

# Each input `.rda` is expected to contain a named object `edgelist`.
# Structure:
# - outer list entries are peptide sequences
# - each sequence contains a list of frames
# - each frame is an edge list for that sequence at one stage / snapshot
# The last frame is treated as the representative assembled network for that sequence/state.
extract_last_frame <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || length(x) == 0) return(NULL)
  frame <- x[[length(x)]]
  if (is.null(frame)) return(NULL)
  frame
}

# Normalize all edge inputs to a 2-column matrix so downstream network construction is consistent.
# Checks here are minimal:
# - allow either data.frame or matrix input
# - require at least 2 columns
# - if extra columns exist, keep only the first 2 as node ids
normalize_edgelist <- function(edge) {
  if (is.null(edge)) return(NULL)
  if (is.data.frame(edge)) edge <- as.matrix(edge)
  if (!is.matrix(edge)) return(NULL)
  if (ncol(edge) < 2) return(NULL)
  if (ncol(edge) > 2) edge <- edge[, 1:2, drop = FALSE]
  edge
}

compute_fractions <- function(edge) {
  edge2 <- normalize_edgelist(edge)
  if (is.null(edge2)) return(NULL)
  # Contacts are modeled as an undirected graph; directionality is not used here.
  net <- network(edge2, directed = FALSE)
  # `fibril_assay()` is bundled locally in `scripts/fibril_assay.R`.
  # It applies the ORCA-based topology classification used in the earlier amyloid-network work
  # and returns node-level membership across ribbon/prism topology classes.
  membership <- fibril_assay(net)
  # Sequence-level topology output is the mean node membership in each class.
  colMeans(membership)
}

# Load the serialized edgelist object and enforce the expected object name and type.
# This validates:
# - file exists
# - the `.rda` actually defines an object named `edgelist`
# - `edgelist` is a non-empty list
load_edgelist_object <- function(rda_path) {
  if (!file.exists(rda_path)) stop('Missing input file: ', rda_path)
  env <- new.env()
  load(rda_path, envir = env)
  if (!exists('edgelist', envir = env)) {
    stop('Expected object `edgelist` in: ', rda_path)
  }
  edgelist <- get('edgelist', envir = env)
  if (!is.list(edgelist) || length(edgelist) == 0) {
    stop('`edgelist` is empty or invalid in: ', rda_path)
  }
  edgelist
}

cat('=== Job Label:', label, '===\n')
cat('Loading monomer:', monomer_file, '\n')
monomer <- load_edgelist_object(monomer_file)
cat('Monomer sequences:', length(monomer), '\n')

cat('Loading dimer single-node:', dimer_file, '\n')
dimer <- load_edgelist_object(dimer_file)
cat('Dimer sequences:', length(dimer), '\n')

common <- intersect(names(monomer), names(dimer))
cat('Common sequences:', length(common), '\n')
if (length(common) == 0) stop('No common sequences between monomer and dimer in ', label)

# Smoke test before the full loop: fail early if the first shared monomer/dimer pair
# cannot be converted into topology fractions.
smoke_seq <- common[[1]]
cat('Smoke test sequence:', smoke_seq, '\n')
mon_frame <- extract_last_frame(monomer[[smoke_seq]])
dim_frame <- extract_last_frame(dimer[[smoke_seq]])
stopifnot(!is.null(mon_frame), !is.null(dim_frame))
mon_frac <- compute_fractions(mon_frame)
dim_frac <- compute_fractions(dim_frame)
if (is.null(mon_frac) || is.null(dim_frac)) stop('Smoke test failed for ', smoke_seq)
cat('Smoke test OK. monomer topologies:', length(mon_frac), 'dimer topologies:', length(dim_frac), '\n')

records <- vector('list', length(common) * 2)
idx <- 1L

for (seq_name in common) {
  # The same peptide sequence is scored independently in monomer and dimer assembly states.
  m_frame <- extract_last_frame(monomer[[seq_name]])
  d_frame <- extract_last_frame(dimer[[seq_name]])

  m_frac <- compute_fractions(m_frame)
  if (!is.null(m_frac)) {
    records[[idx]] <- data.frame(
      source_label = label,
      sequence = seq_name,
      state = 'monomer',
      topology = names(m_frac),
      fraction = as.numeric(m_frac),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  d_frac <- compute_fractions(d_frame)
  if (!is.null(d_frac)) {
    records[[idx]] <- data.frame(
      source_label = label,
      sequence = seq_name,
      state = 'dimer',
      topology = names(d_frac),
      fraction = as.numeric(d_frac),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

records <- records[seq_len(idx - 1L)]
if (length(records) == 0) stop('No valid topology records generated for ', label)

# Final output is a long table with columns:
# - source_label: subset/run identifier
# - sequence: peptide sequence string
# - state: monomer or dimer
# - topology: topology class name from fibril_assay
# - fraction: mean node membership in that class
out_df <- bind_rows(records)
out_dir <- dirname(out_csv)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(out_df, out_csv)

cat('Wrote records:', nrow(out_df), 'to', out_csv, '\n')
