# =============================================================================
# 02_fedfast_scoring.R
#
# Purpose: Build reusable helpers for fed-fasted score protein selection and
#          correlation diagnostics.
#
# Background:
#   Once the fed-fasted processing signature is implicated as the likely cause
#   of unexplained sample clustering, we can select fed-fasted sensitive
#   proteins in two ways:
#     1. Threshold proteins by the absolute SomaLogic Time_12_hours log2 fold
#        change.
#     2. Select less-redundant proteins with a PLINK clumping-like algorithm
#        based on protein-protein correlations.
#
# Inputs:
#   - data/adata_meta.rds, or an existing `adata_meta` object in the current
#     environment. If neither exists, scripts/00_prepare_data.r is sourced.
#   - data/PAV_Plasma_wSpecSheet/fed-fasted_plasma_effect-sizes.csv
#
# Outputs:
#   - fed_fasted                         : cleaned fed-fasted protein table
#   - fed_fasted_cor_all                 : protein correlation matrix using all samples
#   - fed_fasted_cor_reference           : protein correlation matrix excluding
#                                          CLUMP_REFERENCE_EXCLUDE_VALUE
#   - fed_fasted_threshold_proteins      : proteins passing LOGFC_THRESHOLD
#   - fed_fasted_independent_reference   : clumped proteins using reference correlation
#   - fed_fasted_independent_all         : clumped proteins using all-sample correlation
#   - score_from_logfc_threshold()       : calculate score from thresholded proteins
#   - score_from_independent_proteins()  : calculate score from clumped proteins
#   - plot_logfc_threshold_correlation() : plot thresholded-protein correlations
#   - plot_independent_correlation()     : plot color-coded clumping comparison
#   - plot_score_gmm()                   : plot score distribution and GMM fit
#
# Dependencies: tidyverse, corrplot, here, mixtools
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(mixtools))

# =============================================================================
# 0. User-configurable paths and settings
# =============================================================================

ADATA_CACHE <- here("data", "adata_meta.rds")
FED_FASTED_FILE <- here(
  "data",
  "PAV_Plasma_wSpecSheet",
  "fed-fasted_plasma_effect-sizes.csv"
)

# Effect-size threshold used by the focused workflow.
LOGFC_THRESHOLD <- 0.25

# Correlation threshold used for the main independent-protein score example.
CLUMP_COR_THRESHOLD <- 0.75

# Use non-Denver samples as the reference correlation set because the
# fed-fasted-like split is most visible in Denver in this example dataset.
CLUMP_REFERENCE_EXCLUDE_COL <- "site"
CLUMP_REFERENCE_EXCLUDE_VALUE <- "Denver"


# =============================================================================
# 1. Load adata_meta
# =============================================================================

if (!exists("adata_meta")) {
  if (file.exists(ADATA_CACHE)) {
    message("Loading cached adata_meta from ", ADATA_CACHE)
    adata_meta <- readRDS(ADATA_CACHE)
  } else {
    message("adata_meta not found - running 00_prepare_data.r to build it ...")
    source(here("scripts", "00_prepare_data.r"), local = FALSE)
  }
}

if (!exists("adata_meta")) {
  stop(
    "adata_meta was not found after attempting to load or build it.\n",
    "Run scripts/00_prepare_data.r first, or check ADATA_CACHE."
  )
}

if (!CLUMP_REFERENCE_EXCLUDE_COL %in% colnames(adata_meta)) {
  stop(
    "Column `", CLUMP_REFERENCE_EXCLUDE_COL, "` was not found in adata_meta.\n",
    "Update CLUMP_REFERENCE_EXCLUDE_COL in scripts/02_fedfast_scoring.r."
  )
}


# =============================================================================
# 2. Load and clean the fed-fasted sensitive protein table
# =============================================================================

if (!file.exists(FED_FASTED_FILE)) {
  stop(
    "Fed-fasted protein CSV not found:\n  ", FED_FASTED_FILE, "\n",
    "Place fed-fasted_plasma_effect-sizes.csv in data/PAV_Plasma_wSpecSheet/."
  )
}

load_fed_fasted <- function(path, adata_colnames) {
  fed_fasted_raw <- read.csv(path)

  if (!"AptName" %in% colnames(fed_fasted_raw)) {
    fed_fasted_raw <- fed_fasted_raw %>%
      mutate(AptName = paste0("seq.", gsub("-", ".", SeqId)))
  }

  fed_fasted_clean <- fed_fasted_raw %>%
    filter(!is.na(Time_12_hours)) %>%
    filter(AptName %in% adata_colnames) %>%
    mutate(
      abs_Time_12_hours = abs(Time_12_hours),
      Direction = ifelse(Time_12_hours > 0, "Higher after fasting", "Lower after fasting")
    ) %>%
    arrange(desc(abs_Time_12_hours))

  if (nrow(fed_fasted_clean) == 0) {
    stop(
      "No fed-fasted proteins with non-missing Time_12_hours matched adata_meta columns.\n",
      "Check SeqId/AptName formatting in the fed-fasted CSV."
    )
  }

  fed_fasted_clean
}

fed_fasted <- load_fed_fasted(FED_FASTED_FILE, colnames(adata_meta))


# =============================================================================
# 3. Protein selection helpers
# =============================================================================

select_by_logfc_threshold <- function(fed_fasted_df, threshold = LOGFC_THRESHOLD) {
  fed_fasted_df %>%
    filter(abs_Time_12_hours >= threshold)
}

calc_protein_cor <- function(data, proteins) {
  proteins <- intersect(proteins, colnames(data))

  if (length(proteins) < 2) {
    stop("At least two proteins are required to calculate a correlation matrix.")
  }

  data %>%
    select(all_of(proteins)) %>%
    as.data.frame() %>%
    stats::cor(use = "pairwise.complete.obs")
}

get_reference_samples <- function(data) {
  data %>%
    filter(.data[[CLUMP_REFERENCE_EXCLUDE_COL]] != CLUMP_REFERENCE_EXCLUDE_VALUE)
}

get_independent_proteins <- function(fed_fasted_df,
                                     cor_data,
                                     cor_threshold = CLUMP_COR_THRESHOLD) {
  remaining <- fed_fasted_df %>%
    filter(AptName %in% colnames(cor_data)) %>%
    arrange(desc(abs_Time_12_hours))

  independent_proteins <- character()

  while (nrow(remaining) > 0) {
    top_protein <- remaining$AptName[1]
    independent_proteins <- c(independent_proteins, top_protein)

    cor_values <- cor_data[, top_protein]
    keep_proteins <- names(cor_values)[!is.na(cor_values) & cor_values < cor_threshold]

    remaining <- remaining %>%
      filter(AptName %in% keep_proteins) %>%
      arrange(desc(abs_Time_12_hours))
  }

  independent_proteins
}


# =============================================================================
# 4. Correlation matrices and selected protein sets
# =============================================================================

fed_fasted_cor_all <- calc_protein_cor(adata_meta, fed_fasted$AptName)
fed_fasted_cor_reference <- calc_protein_cor(
  get_reference_samples(adata_meta),
  fed_fasted$AptName
)

fed_fasted_threshold <- select_by_logfc_threshold(fed_fasted, LOGFC_THRESHOLD)
fed_fasted_threshold_proteins <- fed_fasted_threshold$AptName

fed_fasted_independent_reference <- get_independent_proteins(
  fed_fasted,
  fed_fasted_cor_reference,
  CLUMP_COR_THRESHOLD
)

fed_fasted_independent_all <- get_independent_proteins(
  fed_fasted,
  fed_fasted_cor_all,
  CLUMP_COR_THRESHOLD
)


# =============================================================================
# 5. Plotting helpers
# =============================================================================

make_combined_cor_matrix <- function(cor_upper,
                                     cor_lower,
                                     fed_fasted_df,
                                     proteins = NULL,
                                     upper_threshold = NULL) {
  if (is.null(proteins)) {
    proteins <- fed_fasted_df$AptName
  }

  order_vec <- fed_fasted_df %>%
    filter(AptName %in% proteins) %>%
    filter(AptName %in% colnames(cor_upper), AptName %in% colnames(cor_lower)) %>%
    arrange(desc(abs_Time_12_hours)) %>%
    pull(AptName)

  if (length(order_vec) < 2) {
    stop("At least two matching proteins are required for a correlation plot.")
  }

  cor_plot <- as.matrix(cor_upper[order_vec, order_vec, drop = FALSE])
  cor_lower_plot <- as.matrix(cor_lower[order_vec, order_vec, drop = FALSE])
  cor_plot[lower.tri(cor_plot)] <- cor_lower_plot[lower.tri(cor_lower_plot)]

  if (!is.null(upper_threshold)) {
    cor_plot[upper.tri(cor_plot) & cor_plot < upper_threshold] <- NA
  }

  fc_lookup <- fed_fasted_df %>%
    filter(AptName %in% order_vec) %>%
    select(AptName, Time_12_hours)

  labels <- vapply(order_vec, function(aptn) {
    fc <- fc_lookup$Time_12_hours[match(aptn, fc_lookup$AptName)]
    paste0(aptn, " (", round(fc, 2), ")")
  }, character(1))

  rownames(cor_plot) <- labels
  colnames(cor_plot) <- labels

  list(matrix = cor_plot, proteins = order_vec)
}

plot_logfc_threshold_correlation <- function(threshold = LOGFC_THRESHOLD,
                                             main = NULL,
                                             label_cex = 0.65) {
  threshold_df <- select_by_logfc_threshold(fed_fasted, threshold)
  plot_data <- make_combined_cor_matrix(
    cor_upper = fed_fasted_cor_all,
    cor_lower = fed_fasted_cor_reference,
    fed_fasted_df = fed_fasted,
    proteins = threshold_df$AptName
  )

  label_cols <- vapply(plot_data$proteins, function(aptn) {
    fc <- fed_fasted$Time_12_hours[match(aptn, fed_fasted$AptName)]
    if (fc > 0) "#B2182B" else "#2166AC"
  }, character(1))

  if (is.null(main)) {
    main <- paste0(
      "Fed-fasted proteins with |Time_12_hours| >= ",
      threshold,
      "\nupper: all samples; lower: excluding ",
      CLUMP_REFERENCE_EXCLUDE_VALUE
    )
  }

  corrplot::corrplot(
    plot_data$matrix,
    tl.col = label_cols,
    tl.cex = label_cex,
    na.label = " ",
    mar = c(0, 0, 3, 0)
  )

  legend(
    "topright",
    legend = c("Higher after fasting", "Lower after fasting"),
    col = c("#B2182B", "#2166AC"),
    pch = 15,
    bty = "n",
    cex = 0.75
  )
  title(main = main, cex.main = 0.95)
}

plot_independent_correlation <- function(cor_threshold = CLUMP_COR_THRESHOLD,
                                         main = NULL,
                                         highlight_cex = 0.7,
                                         other_cex = 0.45) {
  indep_reference <- get_independent_proteins(
    fed_fasted,
    fed_fasted_cor_reference,
    cor_threshold
  )
  indep_all <- get_independent_proteins(
    fed_fasted,
    fed_fasted_cor_all,
    cor_threshold
  )

  plot_data <- make_combined_cor_matrix(
    cor_upper = fed_fasted_cor_all,
    cor_lower = fed_fasted_cor_reference,
    fed_fasted_df = fed_fasted,
    upper_threshold = cor_threshold
  )

  both_sets <- intersect(indep_reference, indep_all)
  only_reference <- setdiff(indep_reference, indep_all)
  only_all <- setdiff(indep_all, indep_reference)

  label_cols <- vapply(plot_data$proteins, function(aptn) {
    if (aptn %in% both_sets) return("purple")
    if (aptn %in% only_reference) return("red")
    if (aptn %in% only_all) return("blue")
    "grey30"
  }, character(1))

  label_cex <- ifelse(
    plot_data$proteins %in% c(indep_reference, indep_all),
    highlight_cex,
    other_cex
  )

  if (is.null(main)) {
    main <- paste0(
      "Independent fed-fasted proteins, cor < ",
      cor_threshold,
      "\nupper: all samples; lower: excluding ",
      CLUMP_REFERENCE_EXCLUDE_VALUE
    )
  }

  corrplot::corrplot(
    plot_data$matrix,
    tl.col = label_cols,
    tl.cex = label_cex,
    na.label = " ",
    mar = c(2.5, 0, 3, 0)
  )

  legend(
    "bottom",
    legend = c(
      paste0("Only excluding ", CLUMP_REFERENCE_EXCLUDE_VALUE, " (", length(only_reference), ")"),
      paste0("Only all samples (", length(only_all), ")"),
      paste0("Both sets (", length(both_sets), ")")
    ),
    col = c("red", "blue", "purple"),
    pch = 15,
    horiz = TRUE,
    xpd = TRUE,
    inset = c(0, -0.02),
    bty = "n",
    cex = 0.75
  )
  title(main = main, cex.main = 0.95)
}


# =============================================================================
# 6. Score and Gaussian mixture model helpers
# =============================================================================

calculate_fedfast_score <- function(fed_fasted_df,
                                    data = adata_meta,
                                    score_name = "Score") {
  score_proteins <- fed_fasted_df %>%
    filter(AptName %in% colnames(data))

  if (nrow(score_proteins) < 2) {
    stop("At least two matched proteins are required to calculate a score.")
  }

  protein_values <- data[, score_proteins$AptName, drop = FALSE] %>%
    as.data.frame()

  centered_values <- sweep(
    protein_values,
    MARGIN = 2,
    STATS = colMeans(protein_values, na.rm = TRUE),
    FUN = "-"
  )

  weighted_values <- sweep(
    centered_values,
    MARGIN = 2,
    STATS = score_proteins$Time_12_hours,
    FUN = "*"
  )

  score <- rowSums(weighted_values, na.rm = TRUE)

  tibble(
    SampleId = data$SampleId,
    !!score_name := score
  )
}

score_from_logfc_threshold <- function(threshold = LOGFC_THRESHOLD,
                                       data = adata_meta,
                                       score_name = "Score") {
  fed_fasted_df <- select_by_logfc_threshold(fed_fasted, threshold)
  calculate_fedfast_score(fed_fasted_df, data = data, score_name = score_name)
}

score_from_independent_proteins <- function(cor_threshold = CLUMP_COR_THRESHOLD,
                                           data = adata_meta,
                                           cor_data = fed_fasted_cor_reference,
                                           score_name = "Score") {
  independent_proteins <- get_independent_proteins(
    fed_fasted,
    cor_data,
    cor_threshold
  )

  fed_fasted_df <- fed_fasted %>%
    filter(AptName %in% independent_proteins)

  calculate_fedfast_score(fed_fasted_df, data = data, score_name = score_name)
}

fit_score_gmm <- function(score_df,
                          score_col = "Score",
                          seed = 100,
                          k = 2) {
  score_values <- score_df[[score_col]]

  if (anyNA(score_values)) {
    stop("Score contains missing values; remove or impute them before GMM fitting.")
  }

  set.seed(seed)
  gmm <- NULL
  invisible(
    utils::capture.output(
      gmm <- mixtools::normalmixEM(
        score_values,
        k = k,
        maxit = 1000,
        epsilon = 1e-08,
        verb = FALSE
      )
    )
  )

  component_order <- order(gmm$mu)
  lower_component <- component_order[1]
  higher_component <- component_order[length(component_order)]
  component <- max.col(gmm$posterior)

  score_df %>%
    mutate(
      GMMClass = case_when(
        component == lower_component ~ "Lower score",
        component == higher_component ~ "Higher score",
        TRUE ~ paste0("Component ", component)
      ),
      GMMClass = factor(GMMClass, levels = c("Lower score", "Higher score"))
    ) %>%
    list(gmm_data = ., gmm = gmm)
}

make_gmm_component_df <- function(gmm, score_range, n = 512) {
  x <- seq(score_range[1], score_range[2], length.out = n)
  component_order <- order(gmm$mu)

  map_dfr(seq_along(gmm$mu), function(component) {
    tibble(
      Score = x,
      Density = gmm$lambda[component] * dnorm(x, gmm$mu[component], gmm$sigma[component]),
      Component = ifelse(
        component == component_order[1],
        "Lower score GMM",
        "Higher score GMM"
      )
    )
  })
}

make_cluster_gmm_label <- function(plot_df) {
  if (!"Cluster" %in% colnames(plot_df)) {
    return(NULL)
  }

  cluster_tab <- plot_df %>%
    count(Cluster, GMMClass) %>%
    tidyr::complete(Cluster, GMMClass, fill = list(n = 0)) %>%
    pivot_wider(names_from = GMMClass, values_from = n, values_fill = 0) %>%
    arrange(Cluster)

  header <- sprintf("%-8s %7s %8s", "Cluster", "Lower", "Higher")
  rows <- sprintf(
    "%-8s %7d %8d",
    as.character(cluster_tab$Cluster),
    cluster_tab$`Lower score`,
    cluster_tab$`Higher score`
  )

  paste(c(header, rows), collapse = "\n")
}

plot_score_gmm <- function(score_df,
                           title = "Fed-fasted score",
                           score_col = "Score",
                           bins = 48,
                           seed = 100) {
  gmm_fit <- fit_score_gmm(score_df, score_col = score_col, seed = seed)
  plot_df <- gmm_fit$gmm_data %>%
    left_join(
      adata_meta %>%
        select(any_of(c("SampleId", "Cluster", "site"))),
      by = "SampleId"
    ) %>%
    rename(Score = all_of(score_col))

  component_df <- make_gmm_component_df(
    gmm_fit$gmm,
    range(plot_df$Score, na.rm = TRUE)
  )

  table_label <- make_cluster_gmm_label(plot_df)
  y_max <- max(component_df$Density, na.rm = TRUE)

  p <- ggplot(plot_df, aes(x = Score)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = bins,
      fill = "white",
      color = "grey45",
      linewidth = 0.25
    ) +
    geom_density(aes(color = Cluster), linewidth = 0.8, adjust = 1.1, na.rm = TRUE) +
    geom_line(
      data = component_df,
      aes(x = Score, y = Density, linetype = Component),
      color = "black",
      linewidth = 0.75
    ) +
    theme_classic() +
    labs(
      title = title,
      x = "Fed-fasted score",
      y = "Density",
      color = "K-means cluster",
      linetype = "GMM component"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )

  if (!is.null(table_label)) {
    p <- p +
      annotate(
        "label",
        x = quantile(plot_df$Score, 0.02, na.rm = TRUE),
        y = y_max * 0.95,
        hjust = 0,
        vjust = 1,
        label = table_label,
        family = "mono",
        size = 3.2
      )
  }

  p
}
