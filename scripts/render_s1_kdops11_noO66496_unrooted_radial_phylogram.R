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
  script_path <- normalizePath("scripts/render_s1_kdops11_noO66496_unrooted_radial_phylogram.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

abs_path <- function(...) file.path(repo_root, ...)

parse_cli <- function(args) {
  cfg <- list(
    tree = abs_path("results", "04_phylogeny_asr", "CoreTree_rooted_MFP_KDOPS11_noO66496.treefile"),
    output_stem = "S1_KDOPS11_noO66496_unrooted_radial_phylogram",
    title = "S1 unrooted radial phylogram",
    subtitle = "KDOPS11 outgroup + MFP | O66496 excluded | branch-length geometry audit",
    layout = Sys.getenv("DAH7PS_LAYOUT_METHOD", unset = "ape"),
    width = 12.0,
    height = 9.0,
    write_log = TRUE,
    backup = TRUE
  )

  if (length(args) == 0L) {
    return(cfg)
  }

  for (arg in args) {
    if (startsWith(arg, "--tree=")) {
      cfg$tree <- normalizePath(sub("^--tree=", "", arg), mustWork = TRUE)
    } else if (startsWith(arg, "--output-stem=")) {
      cfg$output_stem <- sub("^--output-stem=", "", arg)
    } else if (startsWith(arg, "--title=")) {
      cfg$title <- sub("^--title=", "", arg)
    } else if (startsWith(arg, "--subtitle=")) {
      cfg$subtitle <- sub("^--subtitle=", "", arg)
    } else if (startsWith(arg, "--layout=")) {
      cfg$layout <- sub("^--layout=", "", arg)
    } else if (startsWith(arg, "--width=")) {
      cfg$width <- as.numeric(sub("^--width=", "", arg))
    } else if (startsWith(arg, "--height=")) {
      cfg$height <- as.numeric(sub("^--height=", "", arg))
    } else if (arg == "--no-log") {
      cfg$write_log <- FALSE
    } else if (arg == "--no-backup") {
      cfg$backup <- FALSE
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  cfg
}

cfg <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!cfg$layout %in% c("ape", "equal_angle", "daylight")) {
  stop("Unsupported layout: ", cfg$layout, ". Use ape, equal_angle, or daylight.", call. = FALSE)
}
if (!is.finite(cfg$width) || cfg$width <= 0 || !is.finite(cfg$height) || cfg$height <= 0) {
  stop("--width and --height must be positive numbers.", call. = FALSE)
}

data.table::setDTthreads(max(1L, parallel::detectCores(logical = TRUE)))

ia_fasta <- abs_path("results", "02_qc", "nr80_Ia.fasta")
ib_fasta <- abs_path("results", "02_qc", "nr80_Ib.fasta")
ii_fasta <- abs_path("results", "02_qc", "nr80_II.fasta")
figures_dir <- abs_path("figures")
audit_path <- abs_path("results", "04_phylogeny_asr", paste0(cfg$output_stem, "_branch_length_audit.tsv"))
log_path <- abs_path("log.md")

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

MIN_SHADED_CLADE_TIPS <- 25L
MAX_SHADED_CLADES_PER_GROUP <- 24L
DAYLIGHT_MAX_COUNT <- 5L
BRANCH_LENGTH_ABS_TOL <- 1e-7
BRANCH_LENGTH_REL_TOL <- 1e-5

read_fasta_ids <- function(path) {
  lines <- readLines(path, warn = FALSE)
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
  if (nrow(duplicate_ids) > 0L) {
    examples <- paste(head(duplicate_ids$logical_id, 10), collapse = ", ")
    stop("Subtype membership is ambiguous. Examples: ", examples, call. = FALSE)
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

backup_if_exists <- function(path, enabled = TRUE) {
  if (!enabled || !file.exists(path)) {
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

build_layout_data <- function(tree, layout) {
  fortify_phylo <- get("fortify.phylo", envir = asNamespace("ggtree"))
  as.data.frame(fortify_phylo(
    tree,
    layout = layout,
    branch.length = "branch.length",
    MAX_COUNT = DAYLIGHT_MAX_COUNT
  ))
}

audit_branch_length_geometry <- function(tree, tree_data, layout) {
  if (is.null(tree$edge.length)) {
    return(data.table(
      output_stem = cfg$output_stem,
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
      output_stem = cfg$output_stem,
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

  data.table(
    output_stem = cfg$output_stem,
    layout = layout,
    edge_count = nrow(edge_dt),
    max_abs_error = max_abs,
    max_rel_error = max_rel,
    pearson_r = pearson_r,
    audit_status = status,
    scale_bar_status = if (status == "PASS") "pending" else "skipped_audit_fail",
    note = sprintf("abs_tol=%g; rel_tol=%g", BRANCH_LENGTH_ABS_TOL, BRANCH_LENGTH_REL_TOL)
  )
}

build_tip_annotation <- function(tree_data, membership) {
  tips <- as.data.table(tree_data[tree_data$isTip, , drop = FALSE])
  tips[, logical_id := vapply(label, tree_label_to_logical_id, character(1))]
  tips[, plot_group := vapply(label, classify_tip, character(1), membership = membership)]
  tips
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

  invisible(lapply(tree_dt[isTip == FALSE, node], summarize_node))
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
    node_coords <- coords[node %in% c(current_node, descendant_tips)]
    if (nrow(node_coords) < 3L) {
      next
    }
    hull <- node_coords[grDevices::chull(node_coords$x, node_coords$y)]
    polygon_rows[[length(polygon_rows) + 1L]] <- data.table(
      polygon_id = paste(selected$group[[i]], current_node, sep = "_"),
      group = selected$group[[i]],
      x = hull$x,
      y = hull$y
    )
  }

  if (length(polygon_rows) == 0L) data.table() else rbindlist(polygon_rows)
}

make_point_layer <- function(data, mapping = NULL, color = NULL, size = 0.25, shape = 16,
                             alpha = 0.9, show.legend = FALSE) {
  if (nrow(data) == 0L) {
    return(NULL)
  }
  if (requireNamespace("ggrastr", quietly = TRUE)) {
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

nice_scale_width <- function(target) {
  if (!is.finite(target) || target <= 0) {
    return(NA_real_)
  }
  base <- 10^floor(log10(target))
  candidates <- base * c(1, 2, 5, 10)
  candidates[max(which(candidates <= target))]
}

add_scale_bar <- function(p, tree_data, audit) {
  if (!identical(audit$audit_status[[1]], "PASS")) {
    audit$scale_bar_status <- "skipped_audit_fail"
    return(list(plot = p, audit = audit))
  }

  x_rng <- range(tree_data$x, na.rm = TRUE)
  y_rng <- range(tree_data$y, na.rm = TRUE)
  width <- nice_scale_width(diff(x_rng) * 0.10)
  if (!is.finite(width) || width <= 0) {
    audit$scale_bar_status <- "skipped_invalid_width"
    return(list(plot = p, audit = audit))
  }

  scale_x <- x_rng[[1]] + diff(x_rng) * 0.05
  scale_y <- y_rng[[1]] + diff(y_rng) * 0.05
  scale_dt <- data.table(
    x = scale_x,
    xend = scale_x + width,
    y = scale_y,
    label_x = scale_x + width / 2,
    label_y = scale_y - diff(y_rng) * 0.035,
    label = format(width, trim = TRUE, scientific = FALSE, digits = 3)
  )

  p <- p +
    geom_segment(
      data = scale_dt,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      color = "#222222",
      linewidth = 0.42
    ) +
    geom_text(
      data = scale_dt,
      aes(x = label_x, y = label_y, label = label),
      inherit.aes = FALSE,
      color = "#222222",
      size = 3.5
    )
  audit$scale_bar_status <- "added"
  list(plot = p, audit = audit)
}

build_plot <- function(tree, tree_data, tips, edge_dt, shade_dt, audit) {
  palette <- TYPE_COLORS[c("Ia", "Ib", "II", "KDOPS")]
  recognized <- tips[plot_group %in% c("Ia", "Ib", "II")]
  kdops <- tips[plot_group == "KDOPS"]
  other <- tips[plot_group == "Other"]
  audit_label <- if (identical(audit$audit_status[[1]], "PASS")) "PASS" else "FAIL"
  subtitle <- paste0(cfg$subtitle, " ", audit_label, " | UFBoot support not converged")

  p <- ggplot() +
    {
      if (nrow(shade_dt) > 0L) {
        geom_polygon(
          data = shade_dt,
          aes(x = x, y = y, group = polygon_id, fill = group),
          alpha = 0.16,
          color = NA,
          show.legend = FALSE
        )
      } else {
        NULL
      }
    } +
    geom_segment(
      data = edge_dt,
      aes(x = parent_x, y = parent_y, xend = child_x, yend = child_y, color = support_class),
      linewidth = 0.13,
      lineend = "round",
      alpha = 0.92,
      show.legend = TRUE,
      key_glyph = draw_key_path
    ) +
    make_point_layer(
      other,
      mapping = aes(x = x, y = y),
      color = OTHER_COLOR,
      size = 0.13,
      alpha = 0.45,
      show.legend = FALSE
    ) +
    make_point_layer(
      recognized,
      mapping = aes(x = x, y = y, fill = plot_group),
      color = "white",
      size = 0.42,
      shape = 21,
      alpha = 0.92,
      show.legend = TRUE
    ) +
    make_point_layer(
      kdops,
      mapping = aes(x = x, y = y, fill = plot_group),
      color = "#222222",
      size = 1.05,
      shape = 24,
      alpha = 0.98,
      show.legend = TRUE
    ) +
    scale_color_manual(
      values = SUPPORT_COLORS,
      breaks = names(SUPPORT_COLORS),
      limits = names(SUPPORT_COLORS),
      drop = FALSE,
      name = "Branch support"
    ) +
    scale_fill_manual(
      values = palette,
      breaks = names(palette),
      limits = names(palette),
      labels = TYPE_LABELS[names(palette)],
      drop = FALSE,
      name = "Subtype"
    ) +
    labs(title = cfg$title, subtitle = subtitle) +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = c(21, 21, 21, 24),
          size = c(3, 3, 3, 4),
          alpha = 1
        )
      ),
      color = guide_legend(
        override.aes = list(linewidth = c(0.9, 0.9, 0.9), alpha = 1, fill = NA, shape = NA)
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0.13)) +
    scale_y_continuous(expand = expansion(mult = 0.13)) +
    coord_equal(clip = "off") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "#777777", linewidth = 0.28, linetype = "dotted"),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(size = 15, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 8.4, color = "#444444", hjust = 0, margin = margin(b = 5)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 8.5),
      legend.title = element_text(size = 8.5, face = "bold"),
      plot.margin = margin(10, 18, 10, 18)
    )

  add_scale_bar(p, tree_data, audit)
}

save_plot <- function(plot) {
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(figures_dir, paste0(cfg$output_stem, ".pdf"))
  png_path <- file.path(figures_dir, paste0(cfg$output_stem, ".png"))
  pdf_backup <- backup_if_exists(pdf_path, cfg$backup)
  png_backup <- backup_if_exists(png_path, cfg$backup)

  ggsave(
    filename = pdf_path,
    plot = plot,
    width = cfg$width,
    height = cfg$height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )

  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(
    filename = png_path,
    plot = plot,
    width = cfg$width,
    height = cfg$height,
    units = "in",
    dpi = 600,
    device = png_device,
    bg = "white",
    limitsize = FALSE
  )

  data.table(
    pdf = pdf_path,
    png = png_path,
    pdf_backup = pdf_backup,
    png_backup = png_backup
  )
}

git_commit <- function() {
  out <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(err) NA_character_
  )
  if (length(out) == 0L) NA_character_ else out[[1]]
}

append_log <- function(save_dt, audit) {
  if (!cfg$write_log) {
    return(invisible(NULL))
  }

  paths_for_md5 <- c(cfg$tree, save_dt$pdf, save_dt$png, audit_path)
  md5 <- tools::md5sum(paths_for_md5[file.exists(paths_for_md5)])
  md5_lines <- paste0("- `", names(md5), "`: `", unname(as.character(md5)), "`")
  command_text <- paste(
    "mamba run -n dah7ps_ggtree Rscript",
    "scripts/render_s1_kdops11_noO66496_unrooted_radial_phylogram.R"
  )
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  lines <- c(
    "",
    paste0("### ", timestamp, " - R/ggtree S1 KDOPS11 noO66496 unrooted radial phylogram"),
    "",
    paste0("**Command:** `", command_text, "`"),
    paste0("**Commit:** `", git_commit(), "`"),
    paste0("**Input tree:** `", cfg$tree, "`"),
    paste0("**Outputs:** `", save_dt$pdf, "`, `", save_dt$png, "`"),
    paste0("**Audit:** `", audit_path, "`; status=", audit$audit_status[[1]], "; scale_bar=", audit$scale_bar_status[[1]]),
    "**Note:** IQ-TREE UFBoot support did not converge; branch support colors are retained for visual QC and should be interpreted cautiously.",
    "",
    "**MD5:**",
    md5_lines,
    ""
  )
  cat(paste(lines, collapse = "\n"), file = log_path, append = TRUE)
}

message("[render] reading tree: ", cfg$tree)
tree <- ape::read.tree(cfg$tree)
membership <- build_membership()

message("[render] computing ", cfg$layout, " layout")
tree_data <- build_layout_data(tree, cfg$layout)
tips <- build_tip_annotation(tree_data, membership)
edge_dt <- build_edge_segments(tree, tree_data)
audit <- audit_branch_length_geometry(tree, tree_data, cfg$layout)
if (!identical(audit$audit_status[[1]], "PASS")) {
  stop("Branch-length geometry audit failed; refusing to draw certified phylogram.", call. = FALSE)
}

palette <- TYPE_COLORS[c("Ia", "Ib", "II", "KDOPS")]
shade_dt <- make_clade_shading(tree_data, tips, palette)
plot_result <- build_plot(tree, tree_data, tips, edge_dt, shade_dt, audit)
p <- plot_result$plot
audit <- plot_result$audit

invisible(backup_if_exists(audit_path, cfg$backup))
data.table::fwrite(audit, audit_path, sep = "\t", na = "NA")

save_dt <- save_plot(p)
append_log(save_dt, audit)

counts <- tips[, .N, by = plot_group][order(plot_group)]
message("[render] wrote: ", save_dt$pdf)
message("[render] wrote: ", save_dt$png)
message("[render] audit: ", audit_path)
message("[render] tip counts: ", paste(paste0(counts$plot_group, "=", counts$N), collapse = ", "))
