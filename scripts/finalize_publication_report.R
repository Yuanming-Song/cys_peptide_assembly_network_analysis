#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

script_path <- normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- dirname(sub("^--file=", "", script_path))
project_dir <- normalizePath(file.path(script_dir, '..'), mustWork = TRUE)

report_dir <- Sys.getenv("TOPOLOGY_REPORT_DIR", unset = file.path(project_dir, "report"))
plot_dir <- Sys.getenv("TOPOLOGY_PLOT_DIR", unset = file.path(project_dir, "plots"))

state_counts <- read_csv(file.path(report_dir, "dominant_topology_counts.csv"), show_col_types = FALSE)
transition_counts <- read_csv(file.path(report_dir, "transition_counts.csv"), show_col_types = FALSE)

n_seq <- state_counts %>%
  group_by(state) %>%
  summarize(total = sum(count), .groups = "drop")

stopifnot(all(n_seq$total == n_seq$total[1]))
N <- n_seq$total[1]

monomer_none <- state_counts %>% filter(state == "monomer", dominant == "none") %>% pull(percent)
dimer_none <- state_counts %>% filter(state == "dimer", dominant == "none") %>% pull(percent)
monomer_double <- state_counts %>% filter(state == "monomer", dominant == "double 1,2 2-ribbon") %>% pull(percent)
dimer_double <- state_counts %>% filter(state == "dimer", dominant == "double 1,2 2-ribbon") %>% pull(percent)

stay_none <- transition_counts %>% filter(monomer == "none", dimer == "none") %>% pull(count)
stay_double <- transition_counts %>% filter(monomer == "double 1,2 2-ribbon", dimer == "double 1,2 2-ribbon") %>% pull(count)
none_to_double <- transition_counts %>% filter(monomer == "none", dimer == "double 1,2 2-ribbon") %>% pull(count)
double_to_none <- transition_counts %>% filter(monomer == "double 1,2 2-ribbon", dimer == "none") %>% pull(count)

if (length(stay_none) == 0) stay_none <- 0
if (length(stay_double) == 0) stay_double <- 0
if (length(none_to_double) == 0) none_to_double <- 0
if (length(double_to_none) == 0) double_to_none <- 0

pct <- function(x) sprintf("%.1f", 100 * x / N)

# Create publication-style table 1
state_table <- state_counts %>%
  mutate(percent = sprintf("%.1f", percent)) %>%
  arrange(state, desc(count))

lines_state <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Dominant topology composition in monomer and dimer ensembles (N=1000 sequences per state).}",
  "\\label{tab:dominant-composition}",
  "\\begin{tabular}{llrr}",
  "\\toprule",
  "State & Dominant topology & Count & Percent (\\%) \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(state_table))) {
  row <- state_table[i, ]
  lines_state <- c(lines_state, sprintf("%s & %s & %d & %s \\\\", row$state, row$dominant, row$count, row$percent))
}

lines_state <- c(lines_state, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines_state, file.path(report_dir, "table_dominant_composition.tex"))

# Create publication-style table 2 (top transitions)
trans_top <- transition_counts %>%
  mutate(percent = sprintf("%.1f", 100 * count / N)) %>%
  arrange(desc(count)) %>%
  slice_head(n = 10)

lines_trans <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Top monomer-to-dimer dominant-topology transitions.}",
  "\\label{tab:top-transitions}",
  "\\begin{tabular}{llrr}",
  "\\toprule",
  "Monomer dominant & Dimer dominant & Count & Percent (\\%) \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(trans_top))) {
  row <- trans_top[i, ]
  lines_trans <- c(lines_trans, sprintf("%s & %s & %d & %s \\\\", row$monomer, row$dimer, row$count, row$percent))
}

lines_trans <- c(lines_trans, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines_trans, file.path(report_dir, "table_top_transitions.tex"))

manuscript <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{float}",
  "\\usepackage{natbib}",
  "\\title{Topology-State Analysis of Random Peptide Interaction Networks}",
  "\\author{}",
  "\\date{}",
  "\\begin{document}",
  "\\maketitle",
  "\\begin{abstract}",
  sprintf("We analyzed %d peptide sequences with paired monomer and dimer network representations and quantified graphlet-defined amyloid topology classes at the sequence level. Dominant topology distributions shifted between states, with `none` changing from %s\\%% (monomer) to %s\\%% (dimer) and `double 1,2 2-ribbon` from %s\\%% to %s\\%%. Transition analysis identified stable `none` and stable `double 1,2 2-ribbon` as the most frequent outcomes.", N, monomer_none, dimer_none, monomer_double, dimer_double),
  "\\end{abstract}",
  "\\section{Methods}",
  "Edgelists were loaded from `random_sequences_edgelist.rda` (paired monomer/dimer networks). For each network, node-level topology membership was computed with the ORCA-based fibril assay used in prior Grazioli-Butts amyloid network work. Sequence-level topology fractions were computed as class-wise mean node membership, and the dominant class per sequence was assigned as the maximum-fraction class (or `none` if all class fractions were zero; `mixed` for ties).",
  "\\section{Results}",
  sprintf("Each state included %d sequences. The most common dominant topology in monomer was `none` (%s\\%%) followed by `double 1,2 2-ribbon` (%s\\%%). In dimer, `none` (%s\\%%) and `double 1,2 2-ribbon` (%s\\%%) remained dominant.", N, monomer_none, monomer_double, dimer_none, dimer_double),
  sprintf("The largest transition counts were `none \\rightarrow none` (%d, %s\\%%), `double 1,2 2-ribbon \\rightarrow double 1,2 2-ribbon` (%d, %s\\%%), `double 1,2 2-ribbon \\rightarrow none` (%d, %s\\%%), and `none \\rightarrow double 1,2 2-ribbon` (%d, %s\\%%).", stay_none, pct(stay_none), stay_double, pct(stay_double), double_to_none, pct(double_to_none), none_to_double, pct(none_to_double)),
  "\\begin{figure}[H]",
  "\\centering",
  sprintf("\\includegraphics[width=0.95\\linewidth]{%s}", file.path("..", "plots", "topology_fraction_boxplots.png")),
  "\\caption{Distribution of sequence-level topology fractions by state.}",
  "\\end{figure}",
  "\\begin{figure}[H]",
  "\\centering",
  sprintf("\\includegraphics[width=0.85\\linewidth]{%s}", file.path("..", "plots", "dominant_topology_counts.png")),
  "\\caption{Dominant topology counts in monomer and dimer states.}",
  "\\end{figure}",
  "\\begin{figure}[H]",
  "\\centering",
  sprintf("\\includegraphics[width=0.78\\linewidth]{%s}", file.path("..", "plots", "topology_transition_heatmap.png")),
  "\\caption{Monomer-to-dimer transition matrix of dominant topology classes.}",
  "\\end{figure}",
  "\\input{table_dominant_composition}",
  "\\input{table_top_transitions}",
  "\\section{Reproducibility Files}",
  "Core outputs used by this manuscript: `topology_fractions.csv`, `dominant_topology_counts.csv`, `transition_counts.csv`, and the three PNG plots in `../plots`.",
  "\\bibliographystyle{plainnat}",
  "\\bibliography{references}",
  "\\end{document}"
)

writeLines(manuscript, file.path(report_dir, "manuscript.tex"))
cat("Wrote manuscript and publication tables to", report_dir, "\n")
