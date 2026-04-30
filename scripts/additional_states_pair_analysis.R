#!/usr/bin/env Rscript

# Compute paired monomer-vs-dimer network descriptors for each peptide sequence.
# The output `state` column refers specifically to assembly state: monomer or dimer.
# Package roles:
# - readr/dplyr: table assembly and CSV output
# - igraph: graph construction and graph statistics such as assortativity, components,
#   Louvain communities, and coreness
# - network: kept for consistency with the broader project network tooling
# - orca: graphlet orbit counts used to derive role entropy/polarization

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

# Load one serialized edgelist list from an `.rda` file.
# Expected object structure:
# - object name: `edgelist`
# - type: named list
# - each list element: one peptide sequence
# - each sequence element: list of frames across the simulation / processing pipeline
load_edgelist <- function(path) {
  env <- new.env()
  load(path, envir = env)
  if (!exists('edgelist', envir = env)) stop('No edgelist object in ', path)
  x <- get('edgelist', envir = env)
  if (!is.list(x)) stop('edgelist is not list in ', path)
  x
}

# Reduce each sequence trajectory to the last available contact network.
# The last frame is treated as the final assembled-state contact map for that peptide/state.
# Output is forced to a 2-column matrix where each row is one undirected edge (node_i, node_j).
last_frame <- function(x) {
  if (!is.list(x) || length(x) == 0) return(NULL)
  y <- x[[length(x)]]
  if (is.null(y)) return(NULL)
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.matrix(y) || ncol(y) < 2) return(NULL)
  y[,1:2,drop=FALSE]
}

# Basic structural validation before any graph metrics are attempted.
# This function checks only the minimal requirements needed for downstream graph construction:
# - not NULL
# - matrix type
# - exactly 2 columns (source/target node ids)
# - at least 1 edge row
# It does not yet check whether node ids are contiguous/positive; that is handled later
# inside the ORCA-specific remapping step.
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
  # Canonicalize undirected edges before deduplication so (i,j) and (j,i) are treated as the same edge.
  sorted <- t(apply(edge2, 1, function(e) if (e[1] <= e[2]) e else rev(e)))
  unique(sorted)
}

remap_edge_ids <- function(edge2) {
  # `orca::count5()` requires node ids to be contiguous positive integers 1..N.
  # Real project edgelists can have gaps or arbitrary labels, so we remap them here.
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
  # `igraph::graph_from_edgelist()` builds an undirected graph object from the 2-column edge matrix.
  # `igraph::simplify()` then removes self-loops / duplicated edges if any remain.
  g <- graph_from_edgelist(edge2, directed = FALSE)
  simplify(g)
}

# Core-periphery summary is approximated from k-core structure:
# `core_max_k` is the maximum coreness value and `core_frac` is the fraction of nodes in that top core.
core_periphery_score <- function(g) {
  # `igraph::coreness()` returns, for each node, the largest k such that the node belongs
  # to the graph's k-core. Larger values indicate nodes embedded in denser mutually supported structure.
  cvec <- coreness(g)
  if (length(cvec) == 0) return(c(max_core = NA_real_, core_frac = NA_real_))
  kmax <- max(cvec)
  c(max_core = as.numeric(kmax), core_frac = mean(cvec == kmax))
}

# Graphlet-role summary collapses ORCA orbit counts into two scalars:
# entropy for diversity of orbit usage and polarization for concentration into few roles.
graphlet_polarization <- function(edge2) {
  # ORCA is strict: drop non-finite/invalid ids and require positive indices.
  edge2 <- edge2[complete.cases(edge2), , drop = FALSE]
  edge2 <- edge2[edge2[, 1] > 0 & edge2[, 2] > 0, , drop = FALSE]
  if (nrow(edge2) == 0) return(c(role_entropy = NA_real_, role_polarization = NA_real_))
  storage.mode(edge2) <- "integer"
  epos <- remap_edge_ids(edge2)
  if (nrow(epos) == 0) return(c(role_entropy = NA_real_, role_polarization = NA_real_))
  # `orca::count5()` returns node-by-orbit counts for graphlets up to size 5.
  # Summing over nodes gives total usage of each orbit across the whole network.
  od <- count5(epos)
  if (is.null(dim(od))) return(c(role_entropy = NA_real_, role_polarization = NA_real_))
  orb_tot <- colSums(od)
  s <- sum(orb_tot)
  if (s <= 0) return(c(role_entropy = 0, role_polarization = 1))
  p <- orb_tot / s
  p <- p[p > 0]
  # Shannon entropy of orbit usage: higher means roles are spread across more orbit types.
  H <- -sum(p * log(p))
  # `Hmax` is the maximum possible entropy for this orbit vector length.
  Hmax <- log(length(orb_tot))
  # Polarization is reported as 1 - normalized entropy, so higher means more concentrated role usage.
  pol <- ifelse(Hmax > 0, 1 - H / Hmax, NA_real_)
  c(role_entropy = as.numeric(H), role_polarization = as.numeric(pol))
}

# Continuous metrics reported for each monomer or dimer assembly network:
# - assortativity: degree-degree mixing
# - modularity / n_communities: partition structure from Louvain clustering
# - core_max_k / core_frac: strength and size of the highest k-core
# - role_entropy / role_polarization: orbit-role diversity and concentration
# - giant_component_frac: fraction of nodes in the largest connected component
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
  # `igraph::assortativity_degree()` measures degree-degree mixing:
  # positive values mean similar-degree nodes tend to connect;
  # negative values mean high-degree nodes tend to connect to low-degree nodes.
  assort <- suppressWarnings(assortativity_degree(g, directed = FALSE))
  # `igraph::components()` identifies connected components.
  # We keep only the fraction of nodes in the largest component as a cohesion summary.
  comps <- components(g)
  giant_frac <- if (n > 0) max(comps$csize) / n else NA_real_

  mod <- NA_real_
  ncomm <- NA_real_
  if (m > 0 && n > 2) {
    # `igraph::cluster_louvain()` is the Louvain modularity-optimization algorithm.
    # It returns a partition of nodes into communities.
    cl <- cluster_louvain(g)
    # `igraph::modularity()` scores how strongly edges are concentrated within communities
    # relative to a degree-preserving random null.
    mod <- modularity(cl)
    # Number of distinct Louvain communities in the partition.
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

# Sanity check on one shared sequence before processing the full paired set.
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
  # Metrics are computed independently for the monomer and dimer assembly states of the same sequence.
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

# Operational binary flags derived from the continuous descriptors above.
# These are heuristic cutoffs used to summarize whether a network is
# assortative, modular, core-periphery-like, or role-polarized.
# Specifically:
# - `state2_assortative`: assortativity >= 0.10
# - `state2_disassortative`: assortativity <= -0.10
# - `state3_modular`: modularity >= 0.30 and at least 2 communities
# - `state4_core_periphery`: small top-core fraction (<= 0.40) but strong top-core index (>= 3)
# - `state5_role_polarized`: role polarization >= 0.35
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
