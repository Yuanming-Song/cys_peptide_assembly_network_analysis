#!/usr/bin/env Rscript

# Merge per-pair monomer/dimer topology outputs into summary tables and plots.

required_packages <- c('ggplot2','dplyr','tidyr','readr')
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop('Missing packages: ', paste(missing, collapse=', '))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Resolve output locations relative to the repository, but allow overrides for batch execution.
script_path <- normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- dirname(sub("^--file=", "", script_path))
project_dir <- normalizePath(file.path(script_dir, '..'), mustWork = TRUE)

chunk_dir <- Sys.getenv('TOPOLOGY_CHUNK_DIR', unset = file.path(project_dir, 'report', 'chunks'))
report_dir <- Sys.getenv('TOPOLOGY_REPORT_DIR', unset = file.path(project_dir, 'report'))
plot_dir <- Sys.getenv('TOPOLOGY_PLOT_DIR', unset = file.path(project_dir, 'plots'))
run_tag <- Sys.getenv('RUN_TAG', unset = '')
if (run_tag == '') stop('RUN_TAG must be set')

dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

# Each chunk file contains long-format topology fractions for one input pair label.
files <- list.files(chunk_dir, pattern=paste0('^topology_fractions_', run_tag, '_.*\\.csv$'), full.names=TRUE)
if (length(files) == 0) stop('No pair output files found in ', chunk_dir)

fractions_df <- bind_rows(lapply(files, read_csv, show_col_types = FALSE))

# Dominant topology is the maximum-fraction class for a sequence within one assembly state.
label_dominant <- function(df) {
  max_val <- max(df$fraction, na.rm = TRUE)
  if (is.na(max_val) || max_val <= 0) return('none')
  winners <- df$topology[abs(df$fraction - max_val) < 1e-12]
  # A tied maximum is reported as `mixed` rather than forcing a single topology class.
  if (length(winners) == 1) return(winners)
  'mixed'
}

dominant_df <- fractions_df %>%
  group_by(sequence, state) %>%
  summarize(dominant = label_dominant(pick(everything())), .groups='drop')

# Transition table compares the dominant monomer label to the dominant dimer label
# for the same peptide sequence.
transitions <- dominant_df %>%
  pivot_wider(names_from = state, values_from = dominant) %>%
  filter(!is.na(monomer) & !is.na(dimer))

transition_counts <- transitions %>% count(monomer, dimer, name='count') %>% arrange(desc(count))
state_counts <- dominant_df %>% count(state, dominant, name='count') %>% group_by(state) %>% mutate(percent = 100*count/sum(count)) %>% ungroup()

# Persist the merged tables so downstream plotting/reporting can reuse them directly.
write_csv(fractions_df, file.path(report_dir, 'topology_fractions.csv'))
write_csv(state_counts, file.path(report_dir, 'dominant_topology_counts.csv'))
write_csv(transition_counts, file.path(report_dir, 'transition_counts.csv'))

# Distribution plot: how strongly each topology class is represented per sequence/state.
p1 <- ggplot(fractions_df, aes(x=topology, y=fraction, fill=state)) +
  geom_boxplot(outlier.alpha=0.2, position=position_dodge(width=0.8)) +
  labs(title='Topology Fraction Distributions', x='Topology', y='Fraction of nodes') +
  theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1))

ggsave(file.path(plot_dir, 'topology_fraction_boxplots.png'), p1, width=10, height=5, dpi=300)

# Count plot: dominant class frequency in monomer vs dimer assemblies.
p2 <- ggplot(state_counts, aes(x=dominant, y=count, fill=state)) +
  geom_col(position=position_dodge(width=0.8)) +
  labs(title='Dominant Topology Counts', x='Dominant topology', y='Sequence count') +
  theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1))

ggsave(file.path(plot_dir, 'dominant_topology_counts.png'), p2, width=8, height=5, dpi=300)

# Heatmap: sequence-level shift from monomer dominant class to dimer dominant class.
p3 <- ggplot(transition_counts, aes(x=monomer, y=dimer, fill=count)) +
  geom_tile(color='white') +
  geom_text(aes(label=count), size=3) +
  scale_fill_gradient(low='#f0f0f0', high='#2166ac') +
  labs(title='Monomer -> Dimer Topology Transitions', x='Monomer dominant', y='Dimer dominant') +
  theme_bw() + theme(axis.text.x = element_text(angle=30, hjust=1))

ggsave(file.path(plot_dir, 'topology_transition_heatmap.png'), p3, width=7, height=6, dpi=300)

cat('Merged', length(files), 'pair outputs\n')
