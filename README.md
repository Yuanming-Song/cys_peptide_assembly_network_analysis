# CYS Peptide Assembly Network Analysis

This repository contains the core analysis scripts for paired monomer/dimer peptide interaction networks from the CYS peptide screening project.

## Analysis Types

### 1. Topology Assay

`scripts/topology_pair_analysis.R` applies the ORCA-based fibril topology assay used in the earlier Grazioli-Butts amyloid network work. It classifies node participation in motifs such as:

- `1-ribbon`
- `2-ribbon`
- `double 1,2 2-ribbon`
- related prism/ribbon classes
- `none` / `mixed` at the sequence-summary level

For each peptide sequence and state, the script computes topology fractions across nodes. `scripts/topology_merge_from_pairs.R` then merges per-pair outputs into sequence-level summaries and topology plots.

### 2. Additional Network-State Analysis

`scripts/additional_states_pair_analysis.R` computes complementary network-organization metrics that do not assign explicit structural classes, but instead describe how contacts are organized:

- assortativity
- number of communities
- core max-k
- core fraction
- role entropy
- role polarization
- operational state flags such as core-periphery

`scripts/additional_states_distributions_acs_doublecol.R` turns the combined output into a publication-style comparison figure.

## Included Scripts

- `scripts/topology_pair_analysis.R`
- `scripts/topology_merge_from_pairs.R`
- `scripts/additional_states_pair_analysis.R`
- `scripts/additional_states_distributions_acs_doublecol.R`
- `scripts/build_tet_manifest.R`
- `scripts/fibril_assay.R`

The repo does not include cluster submission scripts or LaTeX report sources.

## Input Data

Input edgelists are not bundled here. The expected source is the main project `Edgelist/` tree.

Set:

```bash
export PEPTIDE_DATA_ROOT=/path/to/Edgelist
```

before running manifest generation.

## R Dependencies

- `dplyr`
- `readr`
- `tidyr`
- `ggplot2`
- `igraph`
- `network`
- `orca`

## Minimal Workflow

From the repository root:

1. Build the tetrapeptide manifest.

```bash
Rscript scripts/build_tet_manifest.R
```

2. Run topology analysis for each paired monomer/dimer input.

```bash
Rscript scripts/topology_pair_analysis.R \
  <monomer_rda> <dimer_single_node_rda> <label> \
  report/chunks/topology_fractions_<run_tag>_<label>.csv
```

3. Merge topology outputs.

```bash
RUN_TAG=<run_tag> Rscript scripts/topology_merge_from_pairs.R
```

4. Run additional-states analysis for each paired monomer/dimer input.

```bash
Rscript scripts/additional_states_pair_analysis.R \
  <monomer_rda> <dimer_single_node_rda> <label> \
  report/additional_states/additional_states_metrics_<label>.csv
```

5. Combine per-label additional-states CSVs into:

- `report/additional_states/additional_states_metrics_combined.csv`

6. Generate the publication figure for the additional-state analysis.

```bash
Rscript scripts/additional_states_distributions_acs_doublecol.R
```

## Included Results

This repo includes compact result snapshots under `results/`:

- `results/topology/dominant_topology_counts.csv`
- `results/topology/transition_counts.csv`
- `results/topology/topology_fraction_boxplots.png`
- `results/topology/dominant_topology_counts.png`
- `results/topology/topology_transition_heatmap.png`
- `results/additional_states/summary_flags.csv`
- `results/additional_states/additional_states_metric_distributions_ACS_doublecol.png`

These are small enough to keep under version control and are meant as ready-to-view outputs. Large intermediate chunk files and full per-sequence raw outputs are not bundled.

## Current High-Level Findings

- In the topology assay, `none` and `double 1,2 2-ribbon` are the dominant sequence-level classes in both monomer and dimer states.
- Monomer-to-dimer transitions are dominated by stable `none`, stable `double 1,2 2-ribbon`, and exchanges between those two states.
- In the additional-state analysis, dimerization is associated with:
  - a small decrease in assortativity
  - an increase in community count
  - a drop in core fraction
  - an increase in role entropy
  - a decrease in role polarization
  - a modest increase in the core-periphery flag

These additional metrics are best interpreted as network-organization descriptors rather than direct structural motif labels.
