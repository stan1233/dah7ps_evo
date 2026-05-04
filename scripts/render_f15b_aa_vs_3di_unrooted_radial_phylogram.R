#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_packages <- c("ape", "data.table", "ggplot2", "ggtree", "patchwork")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }

  library(ape)
  library(data.table)
  library(ggplot2)
  library(ggtree)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
if (length(script_file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE)
} else {
  script_path <- normalizePath("scripts/render_f15b_aa_vs_3di_unrooted_radial_phylogram.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

abs_path <- function(...) file.path(repo_root, ...)

aa_tree_path <- abs_path("results", "04_phylogeny_asr", "SkeletonTree_AA.treefile")
threedi_tree_path <- abs_path("results", "04_phylogeny_asr", "SkeletonTree_3Di_Q3Di.treefile")
panel_path <- abs_path("results", "03_msa_modules", "panel35_feature_calibration.tsv")
comparison_path <- abs_path("results", "04_phylogeny_asr", "tree_comparison.tsv")
summary_path <- abs_path("results", "04_phylogeny_asr", "aa_vs_3di_unrooted_radial_summary.tsv")
output_pdf <- abs_path("figures", "F15b_aa_vs_3di_unrooted_radial_phylogram.pdf")
output_png <- abs_path("figures", "F15b_aa_vs_3di_unrooted_radial_phylogram.png")

PDB_ANCHOR_CHAIN <- c(
  "PDB-1KFL" = "PDB-1KFL_A",
  "PDB-1RZM" = "PDB-1RZM_A",
  "PDB-3NV8" = "PDB-3NV8_A",
  "PDB-5CKV" = "PDB-5CKV_A",
  "PDB-2B7O" = "PDB-2B7O_A"
)

PDB_DISPLAY_LABEL <- c(
  "PDB-1KFL_A" = "1KFL*",
  "PDB-1RZM_A" = "1RZM*",
  "PDB-3NV8_A" = "3NV8*",
  "PDB-5CKV_A" = "5CKV*",
  "PDB-2B7O_A" = "2B7O*"
)

SUBTYPE_COLORS <- c(
  "Ia" = "#355C8A",
  "Ib" = "#C46E2B",
  "II" = "#3B8B52"
)

SUBTYPE_LABELS <- c(
  "Ia" = "Type Ia",
  "Ib" = "Type Ib",
  "II" = "Type II"
)

BRANCH_LENGTH_ABS_TOL <- 1e-7
BRANCH_LENGTH_REL_TOL <- 1e-5

backup_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_path <- paste0(path, ".bak_", timestamp)
  ok <- file.rename(path, backup_path)
  if (!ok) {
    stop("Failed to move existing file to backup: ", path, call. = FALSE)
  }
  backup_path
}

extract_accession <- function(rep_id) {
  sub("^UniRef90_", "", rep_id)
}

af_label_from_rep <- function(rep_id) {
  paste0("AF-", extract_accession(rep_id), "-F1-model_v4")
}

tree_label_to_logical_id <- function(label) {
  if (startsWith(label, "AF-") && grepl("-F1-model_v4$", label)) {
    return(paste0("UniRef90_", sub("^AF-(.*)-F1-model_v4$", "\\1", label)))
  }
  if (startsWith(label, "PDB-")) {
    return(strsplit(label, "_", fixed = TRUE)[[1]][[1]])
  }
  label
}

display_label_from_tree_label <- function(label) {
  if (label %in% names(PDB_DISPLAY_LABEL)) {
    return(unname(PDB_DISPLAY_LABEL[[label]]))
  }
  logical_id <- tree_label_to_logical_id(label)
  if (startsWith(logical_id, "UniRef90_")) {
    return(sub("^UniRef90_", "", logical_id))
  }
  logical_id
}

load_panel_mapping <- function() {
  panel <- fread(panel_path, sep = "\t", data.table = FALSE)
  rows <- list()
  for (i in seq_len(nrow(panel))) {
    rep_id <- as.character(panel$rep_id[[i]])
    tree_label <- NA_character_
    if (startsWith(rep_id, "UniRef90_")) {
      tree_label <- af_label_from_rep(rep_id)
    } else if (startsWith(rep_id, "PDB-")) {
      tree_label <- unname(PDB_ANCHOR_CHAIN[[rep_id]])
    }
    if (is.na(tree_label)) {
      next
    }
    rows[[length(rows) + 1L]] <- data.table(
      tree_label = tree_label,
      logical_id = tree_label_to_logical_id(tree_label),
      display_label = display_label_from_tree_label(tree_label),
      subtype = as.character(panel$subtype[[i]]),
      is_anchor = startsWith(rep_id, "PDB-")
    )
  }
  mapping <- rbindlist(rows)
  if (anyDuplicated(mapping$tree_label) > 0) {
    duplicated_labels <- mapping$tree_label[duplicated(mapping$tree_label)]
    stop("Duplicated tree labels in panel mapping: ", paste(duplicated_labels, collapse = ", "), call. = FALSE)
  }
  mapping
}

nice_scale_width <- function(target) {
  if (!is.finite(target) || target <= 0) {
    return(NA_real_)
  }
  base <- 10^floor(log10(target))
  candidates <- base * c(1, 2, 5, 10)
  candidates[max(which(candidates <= target))]
}

set_equal <- function(left, right) {
  length(setdiff(left, right)) == 0L && length(setdiff(right, left)) == 0L
}

node_tip_sets <- function(tree) {
  children <- split(tree$edge[, 2], tree$edge[, 1])
  memo <- new.env(parent = emptyenv())

  collect <- function(node) {
    key <- as.character(node)
    cached <- memo[[key]]
    if (!is.null(cached)) {
      return(cached)
    }
    if (node <= length(tree$tip.label)) {
      result <- tree$tip.label[[node]]
    } else {
      child_nodes <- children[[key]]
      result <- unlist(lapply(child_nodes, collect), use.names = FALSE)
    }
    memo[[key]] <- result
    result
  }

  collect
}

is_unrooted_monophyletic <- function(tree, target_tips) {
  target_tips <- unique(target_tips)
  if (length(target_tips) < 2L) {
    return(NA)
  }
  all_tips <- tree$tip.label
  collect <- node_tip_sets(tree)
  for (child in tree$edge[, 2]) {
    side <- collect(child)
    other <- setdiff(all_tips, side)
    if (set_equal(side, target_tips) || set_equal(other, target_tips)) {
      return(TRUE)
    }
  }
  FALSE
}

maximal_pure_splits <- function(tree, target_tips) {
  target_tips <- unique(target_tips)
  all_tips <- tree$tip.label
  collect <- node_tip_sets(tree)
  pure_sets <- list()

  for (child in tree$edge[, 2]) {
    side <- sort(collect(child))
    other <- sort(setdiff(all_tips, side))
    for (candidate in list(side, other)) {
      if (length(candidate) >= 2L && all(candidate %in% target_tips)) {
        key <- paste(candidate, collapse = "|")
        pure_sets[[key]] <- candidate
      }
    }
  }
  if (length(pure_sets) == 0L) {
    return(list())
  }

  pure_sets <- unname(pure_sets)
  keep <- rep(TRUE, length(pure_sets))
  for (i in seq_along(pure_sets)) {
    for (j in seq_along(pure_sets)) {
      if (i == j) {
        next
      }
      if (length(pure_sets[[j]]) > length(pure_sets[[i]]) && all(pure_sets[[i]] %in% pure_sets[[j]])) {
        keep[[i]] <- FALSE
        break
      }
    }
  }
  pure_sets[keep]
}

load_and_prune_tree <- function(path, mapping) {
  tree <- read.tree(path)
  missing_tips <- setdiff(mapping$tree_label, tree$tip.label)
  if (length(missing_tips) > 0L) {
    stop(
      "Tree is missing panel tip(s): ",
      paste(head(missing_tips, 12), collapse = ", "),
      if (length(missing_tips) > 12L) " ..." else "",
      call. = FALSE
    )
  }
  tree <- keep.tip(tree, mapping$tree_label)
  unroot(tree)
}

build_layout_data <- function(tree) {
  fortify_phylo <- get("fortify.phylo", envir = asNamespace("ggtree"))
  as.data.table(fortify_phylo(
    tree,
    layout = "ape",
    branch.length = "branch.length",
    MAX_COUNT = 5L
  ))
}

build_edge_segments <- function(tree, tree_data) {
  coords <- tree_data[, .(node, x, y)]
  edge_dt <- data.table(
    parent = tree$edge[, 1],
    child = tree$edge[, 2],
    branch_length = if (is.null(tree$edge.length)) NA_real_ else as.numeric(tree$edge.length)
  )
  edge_dt <- merge(edge_dt, coords, by.x = "parent", by.y = "node", all.x = TRUE)
  setnames(edge_dt, c("x", "y"), c("parent_x", "parent_y"))
  edge_dt <- merge(edge_dt, coords, by.x = "child", by.y = "node", all.x = TRUE)
  setnames(edge_dt, c("x", "y"), c("child_x", "child_y"))
  edge_dt
}

audit_branch_length_geometry <- function(edge_dt, tree_id) {
  if (anyNA(edge_dt[, .(parent_x, parent_y, child_x, child_y, branch_length)])) {
    return(data.table(
      tree_id = tree_id,
      edge_count = nrow(edge_dt),
      max_abs_error = NA_real_,
      max_rel_error = NA_real_,
      pearson_r = NA_real_,
      audit_status = "FAIL",
      note = "missing branch length or plotted coordinate"
    ))
  }
  edge_dt[, plotted_length := sqrt((child_x - parent_x)^2 + (child_y - parent_y)^2)]
  edge_dt[, abs_error := abs(plotted_length - branch_length)]
  edge_dt[, rel_error := abs_error / pmax(abs(branch_length), 1e-12)]

  pearson_r <- if (stats::sd(edge_dt$branch_length) == 0 || stats::sd(edge_dt$plotted_length) == 0) {
    NA_real_
  } else {
    stats::cor(edge_dt$branch_length, edge_dt$plotted_length)
  }
  max_abs <- max(edge_dt$abs_error, na.rm = TRUE)
  max_rel <- max(edge_dt$rel_error, na.rm = TRUE)
  status <- if (max_abs <= BRANCH_LENGTH_ABS_TOL || max_rel <= BRANCH_LENGTH_REL_TOL) "PASS" else "FAIL"

  data.table(
    tree_id = tree_id,
    edge_count = nrow(edge_dt),
    max_abs_error = max_abs,
    max_rel_error = max_rel,
    pearson_r = pearson_r,
    audit_status = status,
    note = sprintf("layout=ape; abs_tol=%g; rel_tol=%g", BRANCH_LENGTH_ABS_TOL, BRANCH_LENGTH_REL_TOL)
  )
}

annotate_tips <- function(tree_data, mapping) {
  tips <- tree_data[isTip == TRUE]
  tips <- merge(tips, mapping, by.x = "label", by.y = "tree_label", all.x = TRUE, sort = FALSE)
  if (anyNA(tips$subtype)) {
    stop("Some plotted tips have no subtype annotation.", call. = FALSE)
  }

  extent <- max(diff(range(tree_data$x)), diff(range(tree_data$y)))
  label_pad <- 0.035 * extent
  tips[, radial_distance := sqrt(x^2 + y^2)]
  tips[, unit_x := fifelse(radial_distance > 0, x / radial_distance, 1)]
  tips[, unit_y := fifelse(radial_distance > 0, y / radial_distance, 0)]
  tips[, label_x := x + label_pad * unit_x]
  tips[, label_y := y + label_pad * unit_y]
  tips[, raw_angle := atan2(unit_y, unit_x) * 180 / pi]
  tips[, text_angle := fifelse(raw_angle > 90 | raw_angle < -90, raw_angle + 180, raw_angle)]
  tips[, text_hjust := fifelse(raw_angle > 90 | raw_angle < -90, 1, 0)]
  tips
}

build_plot <- function(tree_id, title, tree, mapping) {
  tree_data <- build_layout_data(tree)
  edge_dt <- build_edge_segments(tree, tree_data)
  audit <- audit_branch_length_geometry(edge_dt, tree_id)
  if (!identical(audit$audit_status[[1]], "PASS")) {
    stop(tree_id, " branch-length geometry audit failed; refusing to draw phylogram.", call. = FALSE)
  }

  tips <- annotate_tips(tree_data, mapping)
  non_anchor_tips <- tips[is_anchor == FALSE]
  anchor_tips <- tips[is_anchor == TRUE]

  edge_range_x <- range(c(edge_dt$parent_x, edge_dt$child_x, tips$x), na.rm = TRUE)
  edge_range_y <- range(c(edge_dt$parent_y, edge_dt$child_y, tips$y), na.rm = TRUE)
  span <- max(diff(edge_range_x), diff(edge_range_y))
  pad <- 0.08 * span
  scale_len <- nice_scale_width(0.16 * span)
  scale_x <- edge_range_x[[1]] + 0.06 * span
  scale_y <- edge_range_y[[1]] - 0.045 * span
  scale_dt <- data.table(
    x = scale_x,
    xend = scale_x + scale_len,
    y = scale_y,
    yend = scale_y,
    label = paste0(scale_len, " subst./site")
  )

  p <- ggplot() +
    geom_segment(
      data = edge_dt,
      aes(x = parent_x, y = parent_y, xend = child_x, yend = child_y),
      linewidth = 0.34,
      color = "#242424",
      lineend = "round"
    ) +
    geom_point(
      data = non_anchor_tips,
      aes(x = x, y = y, color = subtype),
      size = 2.0,
      alpha = 0.96,
      show.legend = TRUE
    ) +
    geom_point(
      data = anchor_tips,
      aes(x = x, y = y, color = subtype),
      shape = 8,
      size = 3.1,
      stroke = 0.9,
      show.legend = FALSE
    ) +
    geom_segment(
      data = scale_dt,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.45,
      color = "#333333"
    ) +
    geom_text(
      data = scale_dt,
      aes(x = (x + xend) / 2, y = y - 0.022 * span, label = label),
      inherit.aes = FALSE,
      size = 2.25,
      color = "#444444"
    ) +
    scale_color_manual(values = SUBTYPE_COLORS, labels = SUBTYPE_LABELS, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3.2, alpha = 1))) +
    coord_equal(
      xlim = c(edge_range_x[[1]] - pad, edge_range_x[[2]] + pad),
      ylim = c(edge_range_y[[1]] - 0.12 * span, edge_range_y[[2]] + pad),
      clip = "off"
    ) +
    labs(title = title, color = NULL) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 6)),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.key.width = unit(0.42, "cm"),
      plot.margin = margin(12, 30, 22, 30)
    )

  list(plot = p, audit = audit, tree_data = tree_data, edge_dt = edge_dt, tips = tips)
}

summarize_tree <- function(tree_id, tree, mapping, audit) {
  type_ii_tips <- mapping[subtype == "II", tree_label]
  type_ii_anchor_tips <- c("PDB-3NV8_A", "PDB-5CKV_A", "PDB-2B7O_A")
  pure_splits <- maximal_pure_splits(tree, type_ii_tips)
  pure_split_sizes <- sort(vapply(pure_splits, length, integer(1)), decreasing = TRUE)
  pure_split_labels <- if (length(pure_splits) == 0L) {
    "none"
  } else {
    paste(vapply(pure_splits, function(items) {
      display <- mapping[match(items, tree_label), display_label]
      paste(display, collapse = ",")
    }, character(1)), collapse = " | ")
  }

  counts <- mapping[, .N, by = subtype]
  count_lookup <- setNames(counts$N, counts$subtype)

  data.table(
    tree_id = tree_id,
    tips = length(tree$tip.label),
    type_Ia_tips = unname(count_lookup[["Ia"]]),
    type_Ib_tips = unname(count_lookup[["Ib"]]),
    type_II_tips = unname(count_lookup[["II"]]),
    all_type_II_unrooted_monophyletic = is_unrooted_monophyletic(tree, type_ii_tips),
    type_II_PDB_anchor_unrooted_monophyletic = is_unrooted_monophyletic(tree, type_ii_anchor_tips),
    type_II_PDB_anchor_group = "O53512_AROG_MYCTU",
    type_II_independent_PDB_anchor_groups = 1L,
    maximal_pure_type_II_split_sizes = if (length(pure_split_sizes) == 0L) "none" else paste(pure_split_sizes, collapse = ";"),
    maximal_pure_type_II_splits = pure_split_labels,
    branch_length_audit_status = audit$audit_status[[1]],
    max_abs_error = audit$max_abs_error[[1]],
    max_rel_error = audit$max_rel_error[[1]]
  )
}

read_nrf <- function() {
  if (!file.exists(comparison_path)) {
    return(NA_real_)
  }
  comparison <- fread(comparison_path, sep = "\t", data.table = FALSE)
  value <- comparison$value[comparison$metric == "RF_normalized"]
  if (length(value) != 1L) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(value))
}

main <- function() {
  for (path in c(aa_tree_path, threedi_tree_path, panel_path)) {
    if (!file.exists(path)) {
      stop("Missing required input: ", path, call. = FALSE)
    }
  }

  mapping <- load_panel_mapping()
  aa_tree <- load_and_prune_tree(aa_tree_path, mapping)
  threedi_tree <- load_and_prune_tree(threedi_tree_path, mapping)

  aa <- build_plot("AA", "AA tree", aa_tree, mapping)
  threedi <- build_plot("3Di", "3Di tree", threedi_tree, mapping)

  nrf <- read_nrf()
  subtitle <- if (is.na(nrf)) {
    "Unrooted radial phylograms pruned to the same 35 logical structure-panel entries"
  } else {
    sprintf(
      "Unrooted radial phylograms pruned to the same 35 logical structure-panel entries | chain-resolved normalized RF = %.4f",
      nrf
    )
  }

  combined <- aa$plot +
    threedi$plot +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(
      title = "AA vs 3Di panel trees",
      subtitle = paste0(subtitle, " | starred Type II PDB tips are O53512 state variants"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 9.5, color = "#555555"),
        plot.margin = margin(8, 8, 8, 8)
      )
    ) &
    theme(legend.position = "bottom")

  summary <- rbindlist(list(
    summarize_tree("AA", aa_tree, mapping, aa$audit),
    summarize_tree("3Di", threedi_tree, mapping, threedi$audit)
  ))

  backup_if_exists(output_pdf)
  backup_if_exists(output_png)
  backup_if_exists(summary_path)

  fwrite(summary, summary_path, sep = "\t")
  ggsave(output_pdf, combined, width = 13.2, height = 7.4, units = "in", bg = "white")
  ggsave(output_png, combined, width = 13.2, height = 7.4, units = "in", dpi = 300, bg = "white")

  message("Wrote: ", output_pdf)
  message("Wrote: ", output_png)
  message("Wrote: ", summary_path)
  message("AA all Type II monophyletic: ", summary[tree_id == "AA", all_type_II_unrooted_monophyletic])
  message("3Di all Type II monophyletic: ", summary[tree_id == "3Di", all_type_II_unrooted_monophyletic])
}

main()
