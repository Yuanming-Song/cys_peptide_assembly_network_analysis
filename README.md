# CYS Peptide Assembly Network Analysis

This directory is a clean subset of the peptide self-assembly analysis workflow.

## What Is Included

- `scripts/topology_pair_analysis.R`
  Computes per-sequence topology fractions from paired monomer/dimer edgelists using the ORCA-based fibril assay.
- `scripts/fibril_assay.R`
  Local copy of the topology classifier required by `topology_pair_analysis.R`.
- `scripts/topology_merge_from_pairs.R`
  Merges per-pair topology outputs into final CSV summaries and topology figures.
- `scripts/additional_states_pair_analysis.R`
  Computes additional network-state metrics from the same paired monomer/dimer edgelists.
- `scripts/additional_states_distributions_acs_doublecol.R`
  Builds the publication figure for additional-state metric distributions.
- `scripts/build_tet_manifest.R`
  Scans the tetrapeptide edgelist tree and builds the manifest of monomer/dimer subset pairs.
- `scripts/finalize_publication_report.R`
  Writes report tables and refreshes `report/manuscript.tex` from merged topology outputs.
- `report/`
  LaTeX report/manuscript sources and generated table fragments.


## Data Dependency

Input data are not bundled here.

Set:

- `PEPTIDE_DATA_ROOT=/path/to/Edgelist`

before running.

## R Package Requirements

Install these packages in the R library you plan to use:

- `dplyr`
- `readr`
- `tidyr`
- `ggplot2`
- `igraph`
- `network`
- `orca`


## Directory Layout

- `report/chunks/`
  Default location for per-pair topology CSV outputs before merge.
- `plots/`
  Default location for generated figures.
- `plots/additional_states/`
  Default location for additional-states figures.

## Minimal Reproduction Flow

1. Build the tetrapeptide manifest.

```bash
Rscript cys_peptide_assembly_network_analysis/scripts/build_tet_manifest.R
```

2. Run topology analysis once per manifest row.

```bash
Rscript cys_peptide_assembly_network_analysis/scripts/topology_pair_analysis.R \
  <monomer_rda> <dimer_single_node_rda> <label> \
  cys_peptide_assembly_network_analysis/report/chunks/topology_fractions_<run_tag>_<label>.csv
```

3. Merge topology outputs into final CSVs and topology plots.

```bash
RUN_TAG=<run_tag> Rscript cys_peptide_assembly_network_analysis/scripts/topology_merge_from_pairs.R
```

4. Run additional-states analysis once per manifest row.

```bash
Rscript cys_peptide_assembly_network_analysis/scripts/additional_states_pair_analysis.R \
  <monomer_rda> <dimer_single_node_rda> <label> \
  cys_peptide_assembly_network_analysis/report/additional_states/additional_states_metrics_<label>.csv
```

5. Combine additional-states per-pair CSVs into a single file named:

- `cys_peptide_assembly_network_analysis/report/additional_states/additional_states_metrics_combined.csv`

6. Generate the additional-states publication figure.

```bash
Rscript cys_peptide_assembly_network_analysis/scripts/additional_states_distributions_acs_doublecol.R
```

7. Refresh manuscript tables and generated manuscript text from topology outputs.

```bash
Rscript cys_peptide_assembly_network_analysis/scripts/finalize_publication_report.R
```

## Environment Variables

- `PEPTIDE_DATA_ROOT`
  Root of the full `Edgelist/` dataset.
- `RUN_TAG`
  Required by `topology_merge_from_pairs.R` to select chunk outputs.
- `TOPOLOGY_CHUNK_DIR`
  Override default topology chunk input directory.
- `TOPOLOGY_REPORT_DIR`
  Override default report output directory.
- `TOPOLOGY_PLOT_DIR`
  Override default plot output directory.
- `ADDITIONAL_STATES_INPUT`
  Override path to `additional_states_metrics_combined.csv`.
- `ADDITIONAL_STATES_PLOT_DIR`
  Override output directory for the additional-states figure.
- `TET_MANIFEST_OUT`
  Override output path for the tetrapeptide manifest CSV.

## Notes

- `topology_pair_analysis.R` depends on `scripts/fibril_assay.R`, which is bundled here so this directory is self-contained at the code level.
- `report/report.tex` and `report/manuscript.tex` expect figures under `plots/` relative to this directory structure.
- If you later publish this as its own repository, this directory can be the repo root directly.
# cys_peptide_assembly_network_analysis
