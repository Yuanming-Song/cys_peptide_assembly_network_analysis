#!/usr/bin/env Rscript

required_packages <- c('readr','dplyr','igraph','network','orca')
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop('Missing packages: ', paste(missing, collapse=', '))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(igraph)
  library(network)
  library(orca)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop('Usage: Rscript states25_pair_analysis.R <monomer_rda> <dimer_rda_single_node> <label> <out_csv>')

monomer_file <- args[[1]]
dimer_file <- args[[2]]
label <- args[[3]]
out_csv <- args[[4]]

load_edgelist <- function(path) {
  env <- new.env()
  load(path, envir = env)
  if (!exists('edgelist', envir = env)) stop('No edgelist object in ', path)
  x <- get('edgelist', envir = env)
  if (!is.list(x)) stop('edgelist is not list in ', path)
  x
}

last_frame <- function(x) {
  if (!is.list(x) || length(x) == 0) return(NULL)
  y <- x[[length(x)]]
  if (is.null(y)) return(NULL)
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.matrix(y) || ncol(y) < 2) return(NULL)
  y[,1:2,drop=FALSE]
}

validate_edge2 <- function(edge2, seq_name, state_name) {
  if (is.null(edge2)) {
    cat('SKIP_NULL', seq_name, state_name, '\n')
    return(FALSE)
  }
  if (!is.matrix(edge2) || ncol(edge2) != 2) {
    cat('SKIP_INVALID_MATRIX', seq_name, state_name, '\n')
    return(FALSE)
  }
  if (nrow(edge2) == 0) {
    cat('SKIP_EMPTY_EDGESET', seq_name, state_name, '\n')
    return(FALSE)
  }
  TRUE
}

deduplicate_edges <- function(edge2) {
  # sort each edge so (i,j) == (j,i) before unique
  sorted <- t(apply(edge2, 1, function(e) if (e[1] <= e[2]) e else rev(e)))
  unique(sorted)
}

remap_edge_ids <- function(edge2) {
  # ORCA::count5 requires strictly positive contiguous integer node IDs.
  clean <- deduplicate_edges(edge2)
  ids <- sort(unique(as.vector(clean)))
  id_map <- setNames(seq_along(ids), ids)
  remapped <- matrix(
    unname(id_map[as.vector(clean)]),
    ncol = 2,
    byrow = FALSE
  )
  storage.mode(remapped) <- 'integer'
  remapped
}

to_igraph <- function(edge2) {
  g <- graph_from_edgelist(edge2, directed = FALSE)
  simplify(g)
}

core_periphery_score <- function(g) {
  cvec <- coreness(g)
  if (length(cvec) == 0) return(c(max_core = NA_real_, core_frac = NA_real_))
  kmax <- max(cvec)
  c(max_core = as.numeric(kmax), core_frac = mean(cvec == kmax))
}

graphlet_polarization <- function(edge2) {
  # ORCA is strict: drop non-finite/invalid ids and require positive indices.
  edge2 <- edge2[complete.cases(edge2), , drop = FALSE]
  edge2 <- edge2[edge2[, 1] > 0 & edge2[, 2] > 0, , drop = FALSE]
  if (nrow(edge2) == 0) return(c(role_entropy = NA_real_, role_polarization = NA_real_))
  storage.mode(edge2) <- "integer"
  epos <- remap_edge_ids(edge2)
  if (nrow(epos) == 0) return(c(role_entropy = NA_real_, role_polarization = NA_real_))
  od <- count5(epos)
  if (is.null(dim(od))) return(c(role_entropy = NA_real_, role_polarization = NA_real_))
  orb_tot <- colSums(od)
  s <- sum(orb_tot)
  if (s <= 0) return(c(role_entropy = 0, role_polarization = 1))
  p <- orb_tot / s
  p <- p[p > 0]
  H <- -sum(p * log(p))
  Hmax <- log(length(orb_tot))
  pol <- ifelse(Hmax > 0, 1 - H / Hmax, NA_real_)
  c(role_entropy = as.numeric(H), role_polarization = as.numeric(pol))
}

compute_metrics <- function(edge2) {
  g <- to_igraph(edge2)
  n <- vcount(g)
  m <- ecount(g)
  if (n < 2 || m < 1) {
    return(data.frame(
      n_nodes = n,
      n_edges = m,
      assortativity = NA_real_,
      modularity = NA_real_,
      n_communities = NA_real_,
      core_max_k = NA_real_,
      core_frac = NA_real_,
      role_entropy = NA_real_,
      role_polarization = NA_real_,
      giant_component_frac = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  assort <- suppressWarnings(assortativity_degree(g, directed = FALSE))
  comps <- components(g)
  giant_frac <- if (n > 0) max(comps$csize) / n else NA_real_

  mod <- NA_real_
  ncomm <- NA_real_
  if (m > 0 && n > 2) {
    cl <- cluster_louvain(g)
    mod <- modularity(cl)
    ncomm <- length(unique(membership(cl)))
  }

  cp <- core_periphery_score(g)
  gp <- graphlet_polarization(edge2)

  data.frame(
    n_nodes = n,
    n_edges = m,
    assortativity = assort,
    modularity = mod,
    n_communities = ncomm,
    core_max_k = cp[['max_core']],
    core_frac = cp[['core_frac']],
    role_entropy = gp[['role_entropy']],
    role_polarization = gp[['role_polarization']],
    giant_component_frac = giant_frac,
    stringsAsFactors = FALSE
  )
}

monomer <- load_edgelist(monomer_file)
dimer <- load_edgelist(dimer_file)
common <- intersect(names(monomer), names(dimer))
if (length(common) == 0) stop('No common sequences for ', label)

# sanity and smoke test
smoke <- common[[1]]
sm_m <- last_frame(monomer[[smoke]])
sm_d <- last_frame(dimer[[smoke]])
if (!validate_edge2(sm_m, smoke, 'monomer') || !validate_edge2(sm_d, smoke, 'dimer')) {
  stop('Smoke test failed edge validation for ', label)
}
sm_mx <- compute_metrics(sm_m)
sm_dx <- compute_metrics(sm_d)
cat('SMOKE_OK', label, smoke, 'mon_assort=', sm_mx$assortativity, 'dim_assort=', sm_dx$assortativity, '\n')

rows <- vector('list', length(common) * 2)
idx <- 1L
for (s in common) {
  mf <- last_frame(monomer[[s]])
  df <- last_frame(dimer[[s]])

  if (validate_edge2(mf, s, 'monomer')) {
    met <- compute_metrics(mf)
    rows[[idx]] <- cbind(data.frame(source_label=label, sequence=s, state='monomer', stringsAsFactors=FALSE), met)
    idx <- idx + 1L
  }
  if (validate_edge2(df, s, 'dimer')) {
    met <- compute_metrics(df)
    rows[[idx]] <- cbind(data.frame(source_label=label, sequence=s, state='dimer', stringsAsFactors=FALSE), met)
    idx <- idx + 1L
  }
}
rows <- rows[seq_len(idx-1L)]
out <- bind_rows(rows)

# operational state flags (2-5)
out <- out %>% mutate(
  state2_assortative = ifelse(!is.na(assortativity) & assortativity >= 0.10, 1L, 0L),
  state2_disassortative = ifelse(!is.na(assortativity) & assortativity <= -0.10, 1L, 0L),
  state3_modular = ifelse(!is.na(modularity) & modularity >= 0.30 & !is.na(n_communities) & n_communities >= 2, 1L, 0L),
  state4_core_periphery = ifelse(!is.na(core_frac) & core_frac <= 0.40 & !is.na(core_max_k) & core_max_k >= 3, 1L, 0L),
  state5_role_polarized = ifelse(!is.na(role_polarization) & role_polarization >= 0.35, 1L, 0L)
)

dir.create(dirname(out_csv), recursive=TRUE, showWarnings=FALSE)
write_csv(out, out_csv)
cat('WROTE', nrow(out), 'rows to', out_csv, '\n')
