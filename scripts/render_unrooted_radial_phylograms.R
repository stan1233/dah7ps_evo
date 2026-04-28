#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_packages <- c("ape", "ggplot2", "ggtree", "data.table")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }

  library(ape)
  library(ggplot2)
  library(ggtree)
  library(data.table)
})

options(stringsAsFactors = FALSE)

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
if (length(script_file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE)
} else {
  script_path <- normalizePath("scripts/render_unrooted_radial_phylograms.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

rel <- function(...) file.path(...)
abs_path <- function(...) file.path(repo_root, ...)

tree_files <- c(
  S1 = rel("results", "04_phylogeny_asr", "CoreTree_rooted_MFP.treefile"),
  S2 = rel("results", "04_phylogeny_asr", "CoreTree_rooted_LGC20.treefile"),
  S3 = rel("results", "04_phylogeny_asr", "CoreTree_rooted_midpoint_ingroup.treefile"),
  S4a = rel("results", "04_phylogeny_asr", "CoreTree_rooted_S4a_top500.treefile"),
  S4b = rel("results", "04_phylogeny_asr", "CoreTree_rooted_S4b_fullsearch.treefile")
)

scenario_ids <- c(
  S1 = "S1_MFP_KDOPS",
  S2 = "S2_LGC20_KDOPS",
  S3 = "S3_MIDPOINT_INGROUP",
  S4a = "S4A_TOP500_PROXY",
  S4b = "S4B_FULLSEARCH_PROXY"
)

ia_fasta <- rel("results", "02_qc", "nr80_Ia.fasta")
ib_fasta <- rel("results", "02_qc", "nr80_Ib.fasta")
ii_fasta <- rel("results", "02_qc", "nr80_II.fasta")
output_dir <- rel("figures")

root_scenarios_path <- rel("results", "04_phylogeny_asr", "root_scenarios.tsv")
audit_path <- rel("results", "04_phylogeny_asr", "unrooted_radial_branch_length_audit.tsv")
log_path <- rel("log.md")

LAYOUT_METHOD <- Sys.getenv("DAH7PS_LAYOUT_METHOD", unset = "ape")
DAYLIGHT_MAX_COUNT <- 5L
STRICT_BRANCH_LENGTH_AUDIT <- TRUE
SHOW_UNKNOWN_TIPS <- TRUE
USE_RASTER_TIPS <- TRUE
MIN_SHADED_CLADE_TIPS <- 25L
MAX_SHADED_CLADES_PER_GROUP <- 24L
PLOT_WIDTH_IN <- 12.0
PLOT_HEIGHT_IN <- 9.0

available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores) || available_cores < 1L) {
  available_cores <- 1L
}
requested_cores <- suppressWarnings(as.integer(Sys.getenv("DAH7PS_PLOT_CORES", unset = NA_character_)))
CPU_CORES <- if (is.na(requested_cores) || requested_cores < 1L) available_cores else min(requested_cores, available_cores)
PARALLEL_WORKERS <- max(1L, min(length(tree_files), CPU_CORES))
WORKER_DT_THREADS <- max(1L, floor(CPU_CORES / PARALLEL_WORKERS))
data.table::setDTthreads(CPU_CORES)

if (!LAYOUT_METHOD %in% c("ape", "equal_angle", "daylight")) {
  stop("Unsupported LAYOUT_METHOD: ", LAYOUT_METHOD, ". Use ape, equal_angle, or daylight.", call. = FALSE)
}

BRANCH_LENGTH_ABS_TOL <- 1e-7
BRANCH_LENGTH_REL_TOL <- 1e-5

TYPE_COLORS <- c(
  Ia = "#355C8A",
  Ib = "#C46E2B",
  II = "#3B8B52",
  KDOPS = "#B6423A"
)
OTHER_COLOR <- "#D2D2D2"
SUPPORT_COLORS <- c(
  ">95% / terminal" = "#202020",
  "90-95%" = "#D6A25A",
  "<90%" = "#C84E5C"
)
TYPE_LABELS <- c(Ia = "Type Ia", Ib = "Type Ib", II = "Type II", KDOPS = "KDOPS")

read_fasta_ids <- function(path) {
  full_path <- abs_path(path)
  lines <- readLines(full_path, warn = FALSE)
  header_lines <- lines[startsWith(lines, ">")]
  ids <- sub("^>", "", header_lines)
  ids <- sub("\\s.*$", "", ids)
  unique(ids[nzchar(ids)])
}

tree_label_to_logical_id <- function(label) {
  label <- as.character(label)
  if (startsWith(label, "AF-") && grepl("-F1-model_v4$", label)) {
    return(paste0("UniRef90_", sub("^AF-(.*)-F1-model_v4$", "\\1", label)))
  }
  if (startsWith(label, "PDB-")) {
    return(strsplit(label, "_", fixed = TRUE)[[1]][[1]])
  }
  label
}

build_membership <- function() {
  subtype_fastas <- c(Ia = ia_fasta, Ib = ib_fasta, II = ii_fasta)
  rows <- rbindlist(lapply(names(subtype_fastas), function(subtype) {
    data.table(logical_id = read_fasta_ids(subtype_fastas[[subtype]]), subtype = subtype)
  }))

  duplicate_ids <- rows[, .(subtypes = paste(sort(unique(subtype)), collapse = ",")), by = logical_id][
    grepl(",", subtypes)
  ]
  if (nrow(duplicate_ids) > 0) {
    examples <- paste(head(duplicate_ids$logical_id, 10), collapse = ", ")
    stop(
      "Subtype membership is ambiguous: IDs found in multiple subtype FASTA files. Examples: ",
      examples,
      call. = FALSE
    )
  }

  membership <- rows$subtype
  names(membership) <- rows$logical_id
  membership
}

classify_tip <- function(label, membership) {
  logical_id <- tree_label_to_logical_id(label)
  subtype <- membership[logical_id]
  if (!is.na(subtype)) {
    return(unname(subtype))
  }
  if (startsWith(logical_id, "KDOPS_")) {
    return("KDOPS")
  }
  "Other"
}

make_tip_annotation <- function(tree_data, membership, scenario) {
  tips <- as.data.table(tree_data[tree_data$isTip, , drop = FALSE])
  tips[, logical_id := vapply(label, tree_label_to_logical_id, character(1))]
  tips[, raw_subtype := vapply(label, classify_tip, character(1), membership = membership)]
  tips[, plot_group := raw_subtype]
  if (!scenario %in% c("S1", "S2")) {
    tips[raw_subtype == "KDOPS", plot_group := "Other"]
  }
  if (!SHOW_UNKNOWN_TIPS) {
    tips[plot_group == "Other", plot_group := NA_character_]
  }

  count_levels <- c("Ia", "Ib", "II", "KDOPS", "Other")
  counts <- as.list(setNames(rep(0L, length(count_levels)), count_levels))
  observed_counts <- table(factor(tips$plot_group, levels = count_levels), useNA = "no")
  for (level in names(observed_counts)) {
    counts[[level]] <- as.integer(observed_counts[[level]])
  }

  raw_kdops <- sum(tips$raw_subtype == "KDOPS", na.rm = TRUE)
  if (!scenario %in% c("S1", "S2") && raw_kdops > 0) {
    warning(
      scenario,
      " contains ",
      raw_kdops,
      " KDOPS-labeled tip(s); plotting them as Other/unexpected.",
      call. = FALSE
    )
  }

  unknown_count <- sum(tips$raw_subtype == "Other", na.rm = TRUE)
  unknown_fraction <- unknown_count / max(nrow(tips), 1L)
  if (unknown_count > 50L || unknown_fraction > 0.01) {
    warning(
      scenario,
      " has ",
      unknown_count,
      " unclassified tip(s) (",
      sprintf("%.2f%%", 100 * unknown_fraction),
      "); plotting them in light grey.",
      call. = FALSE
    )
  }

  list(annotation = tips, counts = counts, raw_kdops = raw_kdops)
}

get_scenario_metadata <- function(scenario) {
  scenarios <- data.table::fread(abs_path(root_scenarios_path), sep = "\t", data.table = FALSE)
  row <- scenarios[scenarios$scenario_id == scenario_ids[[scenario]], , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Expected exactly one root_scenarios.tsv row for ", scenario, call. = FALSE)
  }
  row[1, , drop = FALSE]
}

get_palette_for_scenario <- function(scenario) {
  if (scenario %in% c("S1", "S2")) {
    TYPE_COLORS[c("Ia", "Ib", "II", "KDOPS")]
  } else {
    TYPE_COLORS[c("Ia", "Ib", "II")]
  }
}

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

nice_scale_width <- function(target) {
  if (!is.finite(target) || target <= 0) {
    return(NA_real_)
  }
  base <- 10^floor(log10(target))
  candidates <- base * c(1, 2, 5, 10)
  candidates[max(which(candidates <= target))]
}

audit_branch_length_geometry <- function(tree, tree_data, scenario, layout) {
  if (is.null(tree$edge.length)) {
    return(data.table(
      scenario = scenario,
      layout = layout,
      edge_count = nrow(tree$edge),
      max_abs_error = NA_real_,
      max_rel_error = NA_real_,
      pearson_r = NA_real_,
      audit_status = "FAIL",
      scale_bar_status = "skipped_audit_fail",
      note = "tree has no branch lengths"
    ))
  }
  if (length(tree$edge.length) != nrow(tree$edge) || anyNA(tree$edge.length)) {
    return(data.table(
      scenario = scenario,
      layout = layout,
      edge_count = nrow(tree$edge),
      max_abs_error = NA_real_,
      max_rel_error = NA_real_,
      pearson_r = NA_real_,
      audit_status = "FAIL",
      scale_bar_status = "skipped_audit_fail",
      note = "tree edge.length is missing or contains NA"
    ))
  }

  coords <- as.data.table(tree_data[, c("node", "x", "y"), drop = FALSE])
  edge_dt <- data.table(
    parent = tree$edge[, 1],
    child = tree$edge[, 2],
    branch_length = as.numeric(tree$edge.length)
  )
  edge_dt <- merge(edge_dt, coords, by.x = "parent", by.y = "node", all.x = TRUE)
  setnames(edge_dt, c("x", "y"), c("parent_x", "parent_y"))
  edge_dt <- merge(edge_dt, coords, by.x = "child", by.y = "node", all.x = TRUE)
  setnames(edge_dt, c("x", "y"), c("child_x", "child_y"))

  if (anyNA(edge_dt[, .(parent_x, parent_y, child_x, child_y)])) {
    return(data.table(
      scenario = scenario,
      layout = layout,
      edge_count = nrow(tree$edge),
      max_abs_error = NA_real_,
      max_rel_error = NA_real_,
      pearson_r = NA_real_,
      audit_status = "FAIL",
      scale_bar_status = "skipped_audit_fail",
      note = "missing parent/child coordinates in ggtree layout data"
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
  note <- sprintf("abs_tol=%g; rel_tol=%g", BRANCH_LENGTH_ABS_TOL, BRANCH_LENGTH_REL_TOL)

  data.table(
    scenario = scenario,
    layout = layout,
    edge_count = nrow(edge_dt),
    max_abs_error = max_abs,
    max_rel_error = max_rel,
    pearson_r = pearson_r,
    audit_status = status,
    scale_bar_status = if (status == "PASS") "pending" else "skipped_audit_fail",
    note = note
  )
}

make_point_layer <- function(data, mapping = NULL, color = NULL, size = 0.25, shape = 16,
                             alpha = 0.9, show.legend = FALSE) {
  if (nrow(data) == 0L) {
    return(NULL)
  }
  if (USE_RASTER_TIPS && requireNamespace("ggrastr", quietly = TRUE)) {
    return(ggrastr::geom_point_rast(
      data = data,
      mapping = mapping,
      color = color,
      size = size,
      shape = shape,
      alpha = alpha,
      show.legend = show.legend,
      raster.dpi = 600
    ))
  }
  ggplot2::geom_point(
    data = data,
    mapping = mapping,
    color = color,
    size = size,
    shape = shape,
    alpha = alpha,
    show.legend = show.legend
  )
}

build_layout_data <- function(tree) {
  fortify_phylo <- get("fortify.phylo", envir = asNamespace("ggtree"))
  # The default "ape" method gives the reference-style unrooted phylogram much
  # faster than daylight on 9k-tip trees. "daylight" remains available via
  # DAH7PS_LAYOUT_METHOD=daylight and still records DAYLIGHT_MAX_COUNT.
  as.data.frame(fortify_phylo(
    tree,
    layout = LAYOUT_METHOD,
    branch.length = "branch.length",
    MAX_COUNT = DAYLIGHT_MAX_COUNT
  ))
}

build_edge_segments <- function(tree, tree_data) {
  coords <- as.data.table(tree_data[, c("node", "label", "isTip", "x", "y"), drop = FALSE])
  edge_dt <- data.table(
    parent = tree$edge[, 1],
    child = tree$edge[, 2],
    branch_length = if (is.null(tree$edge.length)) rep(NA_real_, nrow(tree$edge)) else as.numeric(tree$edge.length)
  )
  edge_dt <- merge(edge_dt, coords, by.x = "parent", by.y = "node", all.x = TRUE)
  setnames(edge_dt, c("x", "y", "label", "isTip"), c("parent_x", "parent_y", "parent_label", "parent_isTip"))
  edge_dt <- merge(edge_dt, coords, by.x = "child", by.y = "node", all.x = TRUE)
  setnames(edge_dt, c("x", "y", "label", "isTip"), c("child_x", "child_y", "child_label", "child_isTip"))
  edge_dt[, support := suppressWarnings(as.numeric(child_label))]
  edge_dt[, support_class := fifelse(
    !is.na(support) & support < 90,
    "<90%",
    fifelse(!is.na(support) & support <= 95, "90-95%", ">95% / terminal")
  )]
  edge_dt[, support_class := factor(support_class, levels = names(SUPPORT_COLORS))]
  edge_dt
}

make_clade_shading <- function(tree_data, tips, palette) {
  tree_dt <- as.data.table(tree_data)
  edge_dt <- tree_dt[node != parent, .(parent, child = node)]
  children <- split(edge_dt$child, edge_dt$parent)
  children_names <- names(children)
  tip_group <- tips$plot_group
  names(tip_group) <- as.character(tips$node)
  recognized_groups <- names(palette)
  memo <- new.env(parent = emptyenv())

  summarize_node <- function(node) {
    key <- as.character(node)
    cached <- memo[[key]]
    if (!is.null(cached)) {
      return(cached)
    }
    if (key %in% names(tip_group)) {
      group <- tip_group[[key]]
      if (is.na(group) || !group %in% recognized_groups) {
        group <- "Mixed"
      }
      result <- list(group = group, n_tips = 1L)
      memo[[key]] <- result
      return(result)
    }

    child_nodes <- children[[key]]
    if (is.null(child_nodes) || length(child_nodes) == 0L) {
      result <- list(group = "Mixed", n_tips = 0L)
      memo[[key]] <- result
      return(result)
    }
    child_summaries <- lapply(child_nodes, summarize_node)
    child_groups <- vapply(child_summaries, `[[`, character(1), "group")
    n_tips <- sum(vapply(child_summaries, `[[`, integer(1), "n_tips"))
    clean_groups <- unique(child_groups[child_groups != "Mixed"])
    group <- if (length(clean_groups) == 1L && all(child_groups == clean_groups[[1]])) clean_groups[[1]] else "Mixed"
    result <- list(group = group, n_tips = n_tips)
    memo[[key]] <- result
    result
  }

  internal_nodes <- tree_dt[isTip == FALSE, node]
  invisible(lapply(internal_nodes, summarize_node))

  summary_dt <- rbindlist(lapply(ls(memo), function(key) {
    item <- memo[[key]]
    data.table(node = as.integer(key), group = item$group, n_tips = item$n_tips)
  }))
  summary_dt <- merge(summary_dt, tree_dt[, .(node, parent, isTip)], by = "node", all.x = TRUE)
  parent_summary <- summary_dt[, .(parent = node, parent_group = group)]
  summary_dt <- merge(summary_dt, parent_summary, by = "parent", all.x = TRUE)
  summary_dt[is.na(parent_group), parent_group := "Mixed"]
  selected <- summary_dt[
    isTip == FALSE &
      group %in% recognized_groups &
      n_tips >= MIN_SHADED_CLADE_TIPS &
      group != parent_group
  ]
  if (nrow(selected) == 0L) {
    return(data.table())
  }

  selected <- selected[order(group, -n_tips)]
  selected <- selected[, head(.SD, MAX_SHADED_CLADES_PER_GROUP), by = group]

  collect_tips <- function(node) {
    key <- as.character(node)
    if (key %in% names(tip_group)) {
      return(as.integer(node))
    }
    child_nodes <- children[[key]]
    if (is.null(child_nodes) || length(child_nodes) == 0L) {
      return(integer())
    }
    unlist(lapply(child_nodes, collect_tips), use.names = FALSE)
  }

  polygon_rows <- list()
  coords <- tree_dt[, .(node, x, y)]
  for (i in seq_len(nrow(selected))) {
    current_node <- selected$node[[i]]
    descendant_tips <- collect_tips(current_node)
    target_nodes <- c(current_node, descendant_tips)
    node_coords <- coords[node %in% target_nodes]
    if (nrow(node_coords) < 3L) {
      next
    }
    hull_idx <- grDevices::chull(node_coords$x, node_coords$y)
    hull <- node_coords[hull_idx]
    polygon_id <- paste(selected$group[[i]], current_node, sep = "_")
    polygon_rows[[length(polygon_rows) + 1L]] <- data.table(
      polygon_id = polygon_id,
      group = selected$group[[i]],
      n_tips = selected$n_tips[[i]],
      x = hull$x,
      y = hull$y
    )
  }

  if (length(polygon_rows) == 0L) data.table() else rbindlist(polygon_rows)
}

build_tree_plot <- function(tree, scenario, metadata, membership) {
  message("[", scenario, "] computing ggtree ", LAYOUT_METHOD, " layout")
  tree_data <- build_layout_data(tree)
  tip_info <- make_tip_annotation(tree_data, membership, scenario)
  tips <- tip_info$annotation
  audit <- audit_branch_length_geometry(tree, tree_data, scenario, LAYOUT_METHOD)

  certified <- identical(audit$audit_status[[1]], "PASS")
  if (!certified && STRICT_BRANCH_LENGTH_AUDIT) {
    return(list(
      plot = NULL,
      audit = audit,
      counts = tip_info$counts,
      raw_kdops = tip_info$raw_kdops,
      tree_data = tree_data
    ))
  }

  palette <- get_palette_for_scenario(scenario)
  edge_dt <- build_edge_segments(tree, tree_data)
  shade_dt <- make_clade_shading(tree_data, tips, palette)
  recognized <- tips[plot_group %in% c("Ia", "Ib", "II")]
  kdops <- tips[plot_group == "KDOPS"]
  other <- tips[plot_group == "Other"]

  plot_title <- paste0(scenario, " unrooted radial phylogram")
  plot_subtitle <- if (certified) {
    paste0(as.character(metadata$label[[1]]), " | branch-length geometry audit PASS")
  } else {
    paste0(
      as.character(metadata$label[[1]]),
      " | unrooted radial overview, not branch-length-certified phylogram"
    )
  }

  p <- ggplot2::ggplot() +
    {
      if (nrow(shade_dt) > 0L) {
        ggplot2::geom_polygon(
          data = shade_dt,
          ggplot2::aes(x = x, y = y, group = polygon_id, fill = group),
          alpha = 0.16,
          color = NA,
          show.legend = FALSE
        )
      } else {
        NULL
      }
    } +
    ggplot2::geom_segment(
      data = edge_dt,
      ggplot2::aes(
        x = parent_x,
        y = parent_y,
        xend = child_x,
        yend = child_y,
        color = support_class
      ),
      linewidth = 0.13,
      lineend = "round",
      alpha = 0.92,
      show.legend = TRUE,
      key_glyph = ggplot2::draw_key_path
    ) +
    make_point_layer(
      other,
      mapping = ggplot2::aes(x = x, y = y),
      color = OTHER_COLOR,
      size = 0.13,
      alpha = 0.45,
      show.legend = FALSE
    ) +
    make_point_layer(
      recognized,
      mapping = ggplot2::aes(x = x, y = y, fill = plot_group),
      color = "white",
      size = 0.42,
      shape = 21,
      alpha = 0.92,
      show.legend = TRUE
    ) +
    make_point_layer(
      kdops,
      mapping = ggplot2::aes(x = x, y = y, fill = plot_group),
      color = "#222222",
      size = 1.05,
      shape = 24,
      alpha = 0.98,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = SUPPORT_COLORS,
      breaks = names(SUPPORT_COLORS),
      limits = names(SUPPORT_COLORS),
      drop = FALSE,
      name = "Branch support"
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      breaks = names(palette),
      limits = names(palette),
      labels = TYPE_LABELS[names(palette)],
      drop = FALSE,
      name = "Subtype"
    ) +
    ggplot2::labs(title = plot_title, subtitle = plot_subtitle) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(
          shape = if ("KDOPS" %in% names(palette)) c(21, 21, 21, 24) else c(21, 21, 21),
          size = if ("KDOPS" %in% names(palette)) c(3, 3, 3, 4) else c(3, 3, 3),
          alpha = 1
        )
      ),
      color = ggplot2::guide_legend(override.aes = list(linewidth = c(0.9, 0.9, 0.9), alpha = 1, fill = NA, shape = NA))
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.13)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.13)) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(fill = NA, color = "#777777", linewidth = 0.28, linetype = "dotted"),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 8.4, color = "#444444", hjust = 0, margin = ggplot2::margin(b = 5)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.text = ggplot2::element_text(size = 8.5),
      legend.title = ggplot2::element_text(size = 8.5, face = "bold"),
      plot.margin = ggplot2::margin(10, 18, 10, 18)
    )

  scale_bar_status <- "skipped_audit_fail"
  if (certified) {
    x_rng <- range(tree_data$x, na.rm = TRUE)
    y_rng <- range(tree_data$y, na.rm = TRUE)
    width <- nice_scale_width(diff(x_rng) * 0.10)
    scale_bar_status <- "skipped_invalid_width"
    if (is.finite(width) && width > 0) {
      scale_x <- x_rng[[1]] + diff(x_rng) * 0.05
      scale_y <- y_rng[[1]] + diff(y_rng) * 0.05
      scale_label <- format(width, trim = TRUE, scientific = FALSE, digits = 3)
      scale_dt <- data.table(
        x = scale_x,
        xend = scale_x + width,
        y = scale_y,
        label_x = scale_x + width / 2,
        label_y = scale_y - diff(y_rng) * 0.035,
        label = scale_label
      )
      p <- p +
        ggplot2::geom_segment(
          data = scale_dt,
          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
          inherit.aes = FALSE,
          color = "#222222",
          linewidth = 0.42
        ) +
        ggplot2::geom_text(
          data = scale_dt,
          ggplot2::aes(x = label_x, y = label_y, label = label),
          inherit.aes = FALSE,
          color = "#222222",
          size = 3.5
        )
      scale_bar_status <- "added"
    }
  }
  audit[["scale_bar_status"]] <- scale_bar_status

  list(
    plot = p,
    audit = audit,
    counts = tip_info$counts,
    raw_kdops = tip_info$raw_kdops,
    tree_data = tree_data
  )
}

save_tree_plot <- function(plot, scenario) {
  dir.create(abs_path(output_dir), showWarnings = FALSE, recursive = TRUE)
  stem <- paste0(scenario, "_unrooted_radial_phylogram")
  pdf_path <- abs_path(output_dir, paste0(stem, ".pdf"))
  png_path <- abs_path(output_dir, paste0(stem, ".png"))
  pdf_backup <- backup_if_exists(pdf_path)
  png_backup <- backup_if_exists(png_path)

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = PLOT_WIDTH_IN,
    height = PLOT_HEIGHT_IN,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )

  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = PLOT_WIDTH_IN,
    height = PLOT_HEIGHT_IN,
    units = "in",
    dpi = 600,
    device = png_device,
    bg = "white",
    limitsize = FALSE
  )

  data.table(
    scenario = scenario,
    pdf = pdf_path,
    png = png_path,
    pdf_backup = pdf_backup,
    png_backup = png_backup
  )
}

write_audit_table <- function(audit_dt) {
  dir.create(dirname(abs_path(audit_path)), showWarnings = FALSE, recursive = TRUE)
  backup_if_exists(abs_path(audit_path))
  data.table::fwrite(audit_dt, abs_path(audit_path), sep = "\t", na = "NA")
}

file_md5_table <- function(paths) {
  paths <- unique(paths[file.exists(paths)])
  if (length(paths) == 0L) {
    return(data.table(path = character(), md5 = character()))
  }
  md5 <- tools::md5sum(paths)
  data.table(path = names(md5), md5 = unname(as.character(md5)))
}

git_commit <- function() {
  out <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(err) NA_character_
  )
  if (length(out) == 0L) {
    return(NA_character_)
  }
  out[[1]]
}

package_version_table <- function() {
  pkgs <- c("R", "ape", "ggtree", "treeio", "tidytree", "ggplot2", "data.table", "ggrastr", "ragg")
  data.table(
    package = pkgs,
    version = vapply(pkgs, function(pkg) {
      if (pkg == "R") {
        return(as.character(getRversion()))
      }
      if (requireNamespace(pkg, quietly = TRUE)) {
        return(as.character(utils::packageVersion(pkg)))
      }
      "not_installed"
    }, character(1))
  )
}

markdown_table <- function(dt) {
  if (nrow(dt) == 0L) {
    return(character())
  }
  cols <- names(dt)
  lines <- c(
    paste0("| ", paste(cols, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
  )
  for (i in seq_len(nrow(dt))) {
    values <- vapply(dt[i, ], function(x) {
      value <- as.character(x[[1]])
      if (is.na(value)) "NA" else gsub("\\|", "/", value)
    }, character(1))
    lines <- c(lines, paste0("| ", paste(values, collapse = " | "), " |"))
  }
  lines
}

append_log <- function(output_dt, audit_dt, counts_dt, input_md5_dt, output_md5_dt, pkg_dt, command_text) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  lines <- c(
    "",
    paste0("### ", timestamp, " - R/ggtree unrooted radial phylogram render"),
    "",
    "**Purpose:** Generate one unrooted radial phylogram-style tree per active Phase 4 root scenario with strict branch-length geometry audit.",
    "",
    "**Command:**",
    "",
    "```bash",
    command_text,
    "```",
    "",
    paste0("**Repository commit:** `", git_commit(), "`"),
    paste0("**Layout:** `", LAYOUT_METHOD, "`; `DAYLIGHT_MAX_COUNT=", DAYLIGHT_MAX_COUNT, "` when daylight is requested; `STRICT_BRANCH_LENGTH_AUDIT=", STRICT_BRANCH_LENGTH_AUDIT, "`."),
    paste0("**CPU use:** available cores `", available_cores, "`; requested cores `", CPU_CORES, "`; parallel scenario workers `", PARALLEL_WORKERS, "`; data.table threads per worker `", WORKER_DT_THREADS, "`."),
    paste0("**Reference-style rendering:** black/support-colored branches, translucent subtype clade hulls, no in-plot group labels, and `coord_equal()`."),
    paste0("**Tip rasterization:** `", USE_RASTER_TIPS && requireNamespace("ggrastr", quietly = TRUE), "`; ggrastr version `", package_version_table()[package == "ggrastr", version], "`."),
    "",
    "**Input MD5:**",
    "",
    markdown_table(input_md5_dt),
    "",
    "**Outputs:**",
    "",
    markdown_table(output_dt[, .(scenario, pdf, png)]),
    "",
    "**Output MD5:**",
    "",
    markdown_table(output_md5_dt),
    "",
    "**Subtype Counts:**",
    "",
    markdown_table(counts_dt),
    "",
    "**Branch-Length Geometry Audit:**",
    "",
    markdown_table(audit_dt),
    "",
    "**R Package Versions:**",
    "",
    markdown_table(pkg_dt),
    ""
  )
  cat(paste(lines, collapse = "\n"), file = abs_path(log_path), append = TRUE)
}

validate_inputs <- function() {
  required <- c(tree_files, ia_fasta, ib_fasta, ii_fasta, root_scenarios_path, log_path)
  missing <- required[!file.exists(abs_path(required))]
  if (length(missing) > 0L) {
    stop("Missing required input file(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

render_scenario_worker <- function(scenario, membership) {
  data.table::setDTthreads(WORKER_DT_THREADS)
  tryCatch({
    metadata <- get_scenario_metadata(scenario)
    tree_path <- abs_path(tree_files[[scenario]])
    message("[", scenario, "] reading ", tree_files[[scenario]])
    tree <- ape::read.tree(tree_path)
    if (is.null(tree$edge.length) && STRICT_BRANCH_LENGTH_AUDIT) {
      stop("tree has no branch lengths; strict audit cannot proceed.", call. = FALSE)
    }

    result <- build_tree_plot(tree, scenario, metadata, membership)
    counts <- result$counts
    count_dt <- data.table(
      scenario = scenario,
      Ia = counts$Ia,
      Ib = counts$Ib,
      II = counts$II,
      KDOPS = counts$KDOPS,
      Other = counts$Other,
      raw_KDOPS = result$raw_kdops
    )

    message(
      "[", scenario, "] counts Ia=", counts$Ia,
      " Ib=", counts$Ib,
      " II=", counts$II,
      " KDOPS=", counts$KDOPS,
      " Other=", counts$Other
    )
    message(
      "[", scenario, "] branch audit ",
      result$audit$audit_status[[1]],
      "; max_abs_error=", signif(result$audit$max_abs_error[[1]], 4),
      "; max_rel_error=", signif(result$audit$max_rel_error[[1]], 4)
    )

    if (identical(result$audit$audit_status[[1]], "FAIL") && STRICT_BRANCH_LENGTH_AUDIT) {
      return(list(
        success = FALSE,
        scenario = scenario,
        audit = result$audit,
        counts = count_dt,
        output = data.table(),
        error = paste0("strict branch-length geometry audit failed; no formal plot saved for ", scenario)
      ))
    }

    output_dt <- save_tree_plot(result$plot, scenario)
    message("[", scenario, "] saved ", output_dt$pdf, " and ", output_dt$png)
    list(
      success = TRUE,
      scenario = scenario,
      audit = result$audit,
      counts = count_dt,
      output = output_dt,
      error = NA_character_
    )
  }, error = function(err) {
    list(
      success = FALSE,
      scenario = scenario,
      audit = data.table(),
      counts = data.table(),
      output = data.table(),
      error = paste0("[", scenario, "] ", conditionMessage(err))
    )
  })
}

main <- function() {
  validate_inputs()

  membership <- build_membership()
  input_paths <- c(abs_path(c(tree_files, ia_fasta, ib_fasta, ii_fasta, root_scenarios_path)), script_path)
  input_md5_dt <- file_md5_table(input_paths)
  pkg_dt <- package_version_table()

  audit_rows <- list()
  count_rows <- list()
  output_rows <- list()

  command_text <- "mamba run -n dah7ps_ggtree Rscript scripts/render_unrooted_radial_phylograms.R"
  scenario_names <- names(tree_files)
  message(
    "Parallel rendering: available_cores=", available_cores,
    "; CPU_CORES=", CPU_CORES,
    "; scenario_workers=", PARALLEL_WORKERS,
    "; worker_data_table_threads=", WORKER_DT_THREADS
  )

  if (PARALLEL_WORKERS > 1L) {
    worker_cluster <- parallel::makeCluster(PARALLEL_WORKERS, type = "PSOCK", outfile = "")
    on.exit(parallel::stopCluster(worker_cluster), add = TRUE)
    parallel::clusterExport(worker_cluster, varlist = ls(envir = globalenv()), envir = globalenv())
    parallel::clusterExport(worker_cluster, varlist = "membership", envir = environment())
    parallel::clusterEvalQ(worker_cluster, {
      suppressPackageStartupMessages({
        library(ape)
        library(ggplot2)
        library(ggtree)
        library(data.table)
      })
      data.table::setDTthreads(WORKER_DT_THREADS)
      NULL
    })
    rendered <- parallel::parLapplyLB(
      worker_cluster,
      scenario_names,
      function(scenario) render_scenario_worker(scenario, membership)
    )
  } else {
    rendered <- lapply(scenario_names, render_scenario_worker, membership = membership)
  }
  names(rendered) <- scenario_names

  audit_rows <- lapply(rendered, `[[`, "audit")
  count_rows <- lapply(rendered, `[[`, "counts")
  output_rows <- lapply(rendered, `[[`, "output")
  audit_dt <- rbindlist(audit_rows, fill = TRUE)
  counts_dt <- rbindlist(count_rows, fill = TRUE)
  output_dt <- rbindlist(output_rows, fill = TRUE)
  if (nrow(audit_dt) > 0L) {
    write_audit_table(audit_dt)
  }

  failed <- rendered[!vapply(rendered, `[[`, logical(1), "success")]
  if (length(failed) > 0L) {
    errors <- vapply(failed, `[[`, character(1), "error")
    stop("Renderer did not complete all scenarios:\n", paste(errors, collapse = "\n"), call. = FALSE)
  }

  output_paths <- c(output_dt$pdf, output_dt$png, abs_path(audit_path))
  output_md5_dt <- file_md5_table(output_paths)
  append_log(output_dt, audit_dt, counts_dt, input_md5_dt, output_md5_dt, pkg_dt, command_text)

  message("Audit table: ", abs_path(audit_path))
  message("Log appended: ", abs_path(log_path))
  invisible(0L)
}

if (identical(environment(), globalenv())) {
  main()
}
