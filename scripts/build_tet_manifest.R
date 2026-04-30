#!/usr/bin/env Rscript

# Build manifest of tetrapeptide consolidated subset files.

script_path <- normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- dirname(sub("^--file=", "", script_path))
project_dir <- normalizePath(file.path(script_dir, '..'), mustWork = TRUE)

default_data_root <- normalizePath(file.path(project_dir, '..', 'Edgelist'), mustWork = FALSE)
base <- Sys.getenv('PEPTIDE_DATA_ROOT', unset = default_data_root)
base <- file.path(base, 'Tetrapeptide', 'SubData')
out <- Sys.getenv('TET_MANIFEST_OUT', unset = file.path(project_dir, 'report', 'tet_manifest.csv'))

mono <- list.files(base, pattern='^Tetrapeptide_edgelist_monomer_C[1-4]_[A-Z]\\.rda$', full.names=TRUE)
if (length(mono) == 0) stop('No tetrapeptide monomer subset files found in SubData')

rows <- list()
for (m in mono) {
  b <- basename(m)
  key <- sub('^Tetrapeptide_edgelist_monomer_(C[1-4]_[A-Z])\\.rda$', '\\1', b)
  d <- file.path(base, paste0('Tetrapeptide_edgelist_dimer_', key, '_single_node.rda'))
  if (!file.exists(d)) next
  rows[[length(rows)+1L]] <- data.frame(
    label = paste0('tetra_', key),
    monomer_file = m,
    dimer_file = d,
    stringsAsFactors = FALSE
  )
}

if (length(rows) == 0) stop('No matching monomer/dimer tetrapeptide subset pairs found')
manifest <- do.call(rbind, rows)
manifest <- manifest[order(manifest$label), , drop=FALSE]
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
write.csv(manifest, out, row.names=FALSE, quote=TRUE)
cat('Wrote', nrow(manifest), 'tetrapeptide subset pairs to', out, '\n')
