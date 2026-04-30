#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

script_path <- normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- dirname(sub("^--file=", "", script_path))
project_dir <- normalizePath(file.path(script_dir, '..'), mustWork = TRUE)

input_csv <- Sys.getenv('ADDITIONAL_STATES_INPUT', unset = file.path(project_dir, 'report', 'additional_states', 'additional_states_metrics_combined.csv'))
out_dir <- Sys.getenv('ADDITIONAL_STATES_PLOT_DIR', unset = file.path(project_dir, 'plots', 'additional_states'))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_csv)) stop('Missing combined metrics: ', input_csv)
data <- read_csv(input_csv, show_col_types = FALSE)

metrics_cont <- c('assortativity', 'n_communities', 'core_max_k', 'core_frac', 'role_entropy', 'role_polarization')
metrics_flag <- c('state4_core_periphery')
metrics_plot <- c(metrics_cont, metrics_flag)

collapsed <- data %>%
  group_by(source_label, sequence, state) %>%
  summarize(across(all_of(c(metrics_cont, metrics_flag)), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')

mon <- collapsed %>% filter(state == 'monomer') %>% select(source_label, sequence, all_of(c(metrics_cont, metrics_flag)))
dim <- collapsed %>% filter(state == 'dimer') %>% select(source_label, sequence, all_of(c(metrics_cont, metrics_flag)))
paired <- inner_join(mon, dim, by = c('source_label', 'sequence'), suffix = c('_mon', '_dim'))

if (nrow(paired) == 0) stop('No paired monomer/dimer rows found')

# ACS double-column sizing
width_in <- 7.0
height_in <- 10.5
dpi <- 600

png_file <- file.path(out_dir, 'additional_states_metric_distributions_ACS_doublecol.png')
png(png_file, width = width_in * dpi, height = height_in * dpi, res = dpi)
par(mfrow = c(length(metrics_plot), 2), mar = c(3.2, 1.2, 1.8, 0.8), oma = c(0.8, 2.0, 0.4, 0.2))
par(cex = 0.9)
par(cex.lab = 1.1)
par(cex.axis = 1.0)
par(cex.main = 0.8)
par(mgp = c(2.0, 0.7, 0))

panel_labels <- letters[1:(length(metrics_plot) * 2)]
panel_idx <- 1L
add_panel_label <- function(idx) {
  mtext(sprintf("(%s)", panel_labels[idx]), side = 3, line = 0.1, adj = 0, font = 2, cex = 1.0)
}

for (m in metrics_plot) {
  v_m <- paired[[paste0(m, '_mon')]]
  v_d <- paired[[paste0(m, '_dim')]]
  delta <- v_d - v_m

  if (m == 'state4_core_periphery') {
    xlab_raw <- 'Core–Periphery Flag'
    xlab_delta <- 'Δ Core–Periphery Flag'
    frac_m <- table(factor(v_m, levels = c(0, 1))) / sum(!is.na(v_m))
    frac_d <- table(factor(v_d, levels = c(0, 1))) / sum(!is.na(v_d))
    barplot(
      rbind(frac_m, frac_d),
      beside = TRUE,
      col = c('#3B6EA8', '#B04A4A'),
      names.arg = c('0', '1'),
      main = '',
      ylab = '',
      xlab = xlab_raw,
      legend.text = c('mon', 'dim'),
      args.legend = list(x = 'topright', bty = 'n', cex = 0.8)
    )
    add_panel_label(panel_idx); panel_idx <- panel_idx + 1L

    frac_delta <- table(factor(delta, levels = c(-1, 0, 1))) / sum(!is.na(delta))
    barplot(
      frac_delta,
      col = '#2E7D32',
      names.arg = c('-1', '0', '+1'),
      main = '',
      ylab = '',
      xlab = xlab_delta
    )
    add_panel_label(panel_idx); panel_idx <- panel_idx + 1L
    next
  }

  if (m == 'core_max_k' || m == 'n_communities') {
    xlab_raw <- if (m == 'n_communities') 'Number Of Communities' else 'Core Max-K'
    xlab_delta <- if (m == 'n_communities') 'Δ Number Of Communities' else 'Δ Core Max-K'
    bins <- sort(unique(c(v_m, v_d)))
    bins <- bins[!is.na(bins)]
    frac_m <- table(factor(v_m, levels = bins)) / sum(!is.na(v_m))
    frac_d <- table(factor(v_d, levels = bins)) / sum(!is.na(v_d))
    op <- par(xaxs = 'i')
    barplot(
      rbind(frac_m, frac_d),
      beside = TRUE,
      col = c('#3B6EA8', '#B04A4A'),
      names.arg = bins,
      main = '',
      ylab = '',
      xlab = xlab_raw,
      legend.text = c('mon', 'dim'),
      args.legend = list(x = 'topright', bty = 'n', cex = 0.8),
      xlim = if (m == 'n_communities') c(0, 30) else NULL
    )
    par(op)
    add_panel_label(panel_idx); panel_idx <- panel_idx + 1L

    dbins <- sort(unique(delta))
    dbins <- dbins[!is.na(dbins)]
    frac_delta <- table(factor(delta, levels = dbins)) / sum(!is.na(delta))
    # panel (d): use autoscale (no manual xlim/xticks)
    bp_delta <- barplot(
      frac_delta,
      col = '#2E7D32',
      names.arg = dbins,
      main = '',
      ylab = '',
      xlab = xlab_delta
    )
    if (m == 'n_communities') {
      panel_d_xlim <- par('usr')[1:2]
      cat('PANEL_D_XLIM', panel_d_xlim[1], panel_d_xlim[2], '\\n')
    }
    add_panel_label(panel_idx); panel_idx <- panel_idx + 1L
    next
  }

  dm <- density(v_m, na.rm = TRUE)
  dd <- density(v_d, na.rm = TRUE)
  ylim <- range(0, dm$y, dd$y)
  if (m == 'role_entropy') {
    plot(dm, main = '', xlab = 'Role Entropy', ylab = '', col = '#3B6EA8', lwd = 1.2, ylim = ylim, xlim = c(0, 5))
  } else if (m == 'role_polarization') {
    plot(dm, main = '', xlab = 'Role Polarization', ylab = '', col = '#3B6EA8', lwd = 1.2, ylim = ylim, xlim = c(0, 1))
  } else if (m == 'n_communities') {
    xr <- range(c(v_m, v_d), na.rm = TRUE)
    plot(dm, main = '', xlab = 'Number Of Communities', ylab = '', col = '#3B6EA8', lwd = 1.2, ylim = ylim, xlim = xr)
    rug(v_m, col = '#3B6EA8', ticksize = 0.02)
    rug(v_d, col = '#B04A4A', ticksize = 0.02)
  } else {
    xr <- range(c(v_m, v_d), na.rm = TRUE)
    xlab_raw <- if (m == 'assortativity') 'Assortativity' else if (m == 'core_frac') 'Core Fraction' else m
    plot(dm, main = '', xlab = xlab_raw, ylab = '', col = '#3B6EA8', lwd = 1.2, ylim = ylim, xlim = xr)
  }
  lines(dd, col = '#B04A4A', lwd = 1.2)
  add_panel_label(panel_idx); panel_idx <- panel_idx + 1L

  ddel <- density(delta, na.rm = TRUE)
  if (m == 'role_entropy') {
    plot(ddel, main = '', xlab = 'Δ Role Entropy', ylab = '', col = '#2E7D32', lwd = 1.2, xlim = c(-5, 5))
  } else if (m == 'role_polarization') {
    plot(ddel, main = '', xlab = 'Δ Role Polarization', ylab = '', col = '#2E7D32', lwd = 1.2, xlim = c(-1, 1))
  } else {
    xr_d <- range(delta, na.rm = TRUE)
    xlab_delta <- if (m == 'assortativity') 'Δ Assortativity' else if (m == 'core_frac') 'Δ Core Fraction' else paste('Δ', m)
    plot(ddel, main = '', xlab = xlab_delta, ylab = '', col = '#2E7D32', lwd = 1.2, xlim = xr_d)
  }
  add_panel_label(panel_idx); panel_idx <- panel_idx + 1L
}

mtext('Density', side = 2, outer = TRUE, line = 0.6, cex = 1.0)

dev.off()
cat('WROTE', png_file, '\\n')
