# =============================================================================
# 01_cluster_diagnosis.R
#
# Purpose: Diagnose the cause of unexplained sample clustering in SomaScan
#          proteomics data by comparing between-cluster differential expression
#          to SomaLogic-provided lists of processing-sensitive proteins.
#
# Background:
#   After PCA and k-means clustering, two sample clusters may emerge that
#   cannot be explained by experimental variables (e.g., disease status, sex,
#   age, ancestry). SomaLogic provides lists of proteins known to be sensitive
#   to specific pre-analytical processing differences:
#     - Fed vs. fasted state
#     - Freeze-thaw cycles
#     - Time to spin
#     - Time to decant
#     - Time to freeze
#
#   If one of these lists shows strong concordance between its known
#   processing fold-changes and the observed between-cluster effect sizes,
#   that processing variable is likely responsible for the clustering.
#
# Inputs:
#   - adata_meta : SomaScan ADAT data joined with sample metadata,
#                  log2-transformed, with a 'Cluster' column ('Big'/'Small')
#                  added from a prior k-means step.
#                  Must be present in the environment before this script is
#                  sourced, OR a cached copy must exist at
#                  data/de_cluster_results.rds (produced on first run).
#
#   - data/      : Directory containing SomaLogic-provided CSVs, one per
#                  processing condition. Files are gitignored and must be
#                  obtained directly from SomaLogic. Expected filenames
#                  (adjust the pattern below if yours differ):
#                    FedFasted.csv, FreezeThaw.csv, TimeToSpin.csv,
#                    TimeToDecant.csv, TimeToFreeze.csv
#                  Each CSV must contain at minimum:
#                    SeqId        - SomaScan aptamer identifier
#                    [effect col] - log2 fold change at the highest available
#                                   time point or cycle count
#                                   (e.g. Time_12_hours, cycle_10, Time_24_hours)
#
# Outputs:
#   - data/de_cluster_results.rds  : cached per-protein linear regression
#                                    results (written on first run, reused
#                                    on subsequent runs)
#   - diag_figure                  : ggplot object (5-panel diagnostic figure)
#                                    returned to the calling environment
#
# User-configurable settings: see Section 0 below.
#
# Dependencies: tidyverse, SomaDataIO, ggpmisc, ggpubr, here
# =============================================================================

library(tidyverse)
library(SomaDataIO)
library(ggpmisc)
library(here)

# =============================================================================
# 0. User-configurable settings
#    Adjust these to match your dataset before running.
# =============================================================================

# Column in adata_meta that identifies the clustered subset of samples.
# Set to NULL to use all samples.
CLUSTER_SUBSET_COL   <- "site"       # e.g. "site", "cohort", NULL
CLUSTER_SUBSET_VALUE <- "Denver"     # value in CLUSTER_SUBSET_COL to keep

# Covariates to include in the linear model alongside Cluster.
# Remove or add column names to match your metadata.
COVARIATES <- c(
  "Active_Asthma",   # primary phenotype — adjust name to yours
  "Sex",
  "Age",
  "full_PC1",        # ancestry PC1
  "full_PC2",        # ancestry PC2
  "SomascanBatch"    # assay batch
)

# Columns to factorize before model fitting
FACTOR_COLS <- c("Cluster", "Sex", "SomascanBatch")

# Path to cached DE results (written on first run, loaded on subsequent runs)
CACHE_PATH <- here("data", "de_cluster_results.rds")

# Number of top positive and negative processing-effect proteins to show in
# each diagnostic panel, based on the x-axis effect-size column. Use Inf to
# plot all proteins.
TOP_PROCESSING_PROTEINS <- 5


# =============================================================================
# 1. Load the processing-sensitive protein lists
# =============================================================================
# SomaLogic provides these CSVs upon request.
# Place them in:  data/PAV_Plasma_wSpecSheet/
# They are excluded from version control via .gitignore.
# Contact SomaLogic to obtain them.
#
# Expected filenames (the pattern below matches these):
#   FedFasted.csv, FreezeThaw.csv, TimeToSpin.csv,
#   TimeToDecant.csv, TimeToFreeze.csv

CSV_DIR <- here("data", "PAV_Plasma_wSpecSheet")

sens_prot_files <- list.files(
  CSV_DIR,
  pattern    = "effect.sizes?\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(sens_prot_files) == 0) {
  stop(
    "No processing-sensitive protein CSVs found in data/PAV_Plasma_wSpecSheet/.\n",
    "Expected files: fed-fasted_plasma_effect-sizes.csv, ",
    "freeze-thaw_plasma_effect-sizes.csv, time-to-decant_plasma_effect_sizes.csv, ",
    "time-to-freeze_plasma_effect-sizes.csv, time-to-spin_plasma_effect_sizes.csv.\n",
    "Obtain these from SomaLogic and place them in data/PAV_Plasma_wSpecSheet/."
  )
}

# Assign short friendly names based on the condition mentioned in each filename
sens_prot_names <- dplyr::case_when(
  grepl("fed.fasted",   basename(sens_prot_files), ignore.case = TRUE) ~ "FedFasted",
  grepl("freeze.thaw",  basename(sens_prot_files), ignore.case = TRUE) ~ "FreezeThaw",
  grepl("time.to.decant", basename(sens_prot_files), ignore.case = TRUE) ~ "TimeToDecant",
  grepl("time.to.freeze", basename(sens_prot_files), ignore.case = TRUE) ~ "TimeToFreeze",
  grepl("time.to.spin",   basename(sens_prot_files), ignore.case = TRUE) ~ "TimeToSpin",
  TRUE ~ tools::file_path_sans_ext(basename(sens_prot_files))
)
names(sens_prot_files) <- sens_prot_names
sens_prot_list <- lapply(sens_prot_files, read.csv)


# =============================================================================
# 2. Fit per-protein linear models (or load cached results)
# =============================================================================
# For each human protein, regress log2 RFU on Cluster assignment while
# adjusting for the covariates defined in Section 0. A positive beta
# coefficient for ClusterSmall means the protein is more abundant in the
# Small cluster.
#
# Model fitting is skipped if a cached RDS exists at CACHE_PATH.

if (file.exists(CACHE_PATH)) {

  message("Loading cached DE results from ", CACHE_PATH)
  de_cluster_results <- readRDS(CACHE_PATH)

} else {

  if (!exists("adata_meta")) {
    message("adata_meta not found — running 00_prepare_data.r to build it ...")
    source(here("scripts", "00_prepare_data.r"), local = FALSE)
  }

  if (!exists("adata_meta")) {
    stop(
      "adata_meta still not found after running 00_prepare_data.r.\n",
      "Check that ADAT_FILE, PHENO_FILE, and SHIPMENT_FILE are correct in\n",
      "scripts/00_prepare_data.r."
    )
  }

  # --- Subset to the clustered samples ---
  if (!is.null(CLUSTER_SUBSET_COL)) {
    adata_subset <- adata_meta %>%
      filter(.data[[CLUSTER_SUBSET_COL]] == CLUSTER_SUBSET_VALUE)
  } else {
    adata_subset <- adata_meta
  }

  # Factorize specified columns
  for (col in intersect(FACTOR_COLS, colnames(adata_subset))) {
    adata_subset[[col]] <- factor(adata_subset[[col]])
  }
  # Ensure 'Big' is the reference level so a positive beta = higher in Small
  if ("Cluster" %in% colnames(adata_subset)) {
    adata_subset$Cluster <- relevel(factor(adata_subset$Cluster), ref = "Big")
  }

  # --- Retrieve human analyte metadata ---
  analytes <- getAnalyteInfo(adat = adata_subset) %>%
    filter(Organism == "Human") %>%
    select(AptName, SeqId, Target, UniProt, EntrezGeneID, Dilution)

  # --- Build formula for each protein ---
  covariate_str <- paste(COVARIATES, collapse = " + ")
  analytes <- analytes %>%
    mutate(
      formula = map(
        AptName,
        ~ as.formula(paste0(.x, " ~ Cluster + ", covariate_str))
      ),
      model = map(formula, ~ lm(.x, data = adata_subset))
    )

  # --- Extract Cluster beta, p-value, and R² from each model ---
  extract_cluster_coef <- function(apt_name, model) {
    coefs <- summary(model)$coefficients
    coefs <- as.data.frame(coefs)
    colnames(coefs) <- c("beta", "se", "t", "p.value")
    coefs <- tibble::rownames_to_column(coefs, "Variable")

    coefs %>%
      filter(str_detect(Variable, "Cluster")) %>%
      select(Variable, beta, p.value) %>%
      mutate(
        AptName = apt_name,
        R2      = summary(model)$r.squared
      )
  }

  model_params <- map2_dfr(
    analytes$AptName,
    analytes$model,
    extract_cluster_coef
  ) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

  de_cluster_results <- analytes %>%
    select(AptName, SeqId, Target, UniProt, EntrezGeneID, Dilution) %>%
    left_join(model_params, by = "AptName") %>%
    mutate(direction = ifelse(beta > 0, "HighInSmall", "LowInSmall"))

  # Cache results so this step can be skipped on re-render
  saveRDS(de_cluster_results, CACHE_PATH)
  message("DE results saved to ", CACHE_PATH)
}


# =============================================================================
# 3. Join processing-sensitive protein lists with cluster DE results
# =============================================================================
# Merge the known processing effect size (x-axis) with the observed
# between-cluster beta (y-axis). Concordance in quadrants I (+/+) and
# III (-/-) of the resulting plot indicates that processing condition is
# driving the sample clustering.

sens_prot_list <- lapply(sens_prot_list, function(df) {
  df %>% left_join(de_cluster_results, by = "SeqId")
})


# =============================================================================
# 4. Build the 5-panel diagnostic figure
# =============================================================================
# Each panel corresponds to one processing condition. The x-axis is the
# SomaLogic-provided log2 fold change; the y-axis is the observed cluster
# beta. Point size and colour encode -log10 adjusted p-value. Red crosshairs
# mark the origin. Concordance in quadrants I and III is the diagnostic signal.

make_diag_panel <- function(data, x_col, title_label, top_processing_proteins = Inf) {
  plot_all_processing_proteins <- is.infinite(top_processing_proteins) &&
    top_processing_proteins > 0

  if (
    !is.numeric(top_processing_proteins) ||
      length(top_processing_proteins) != 1 ||
      is.na(top_processing_proteins) ||
      (!plot_all_processing_proteins &&
         (top_processing_proteins < 1 || top_processing_proteins %% 1 != 0))
  ) {
    stop("top_processing_proteins must be a positive whole number or Inf.")
  }

  # Drop rows where either axis value is missing
  data <- data %>% filter(!is.na(.data[[x_col]]), !is.na(beta))

  if (!plot_all_processing_proteins) {
    data <- bind_rows(
      data %>%
        filter(.data[[x_col]] > 0) %>%
        slice_max(
          order_by = .data[[x_col]],
          n = top_processing_proteins,
          with_ties = FALSE
        ),
      data %>%
        filter(.data[[x_col]] < 0) %>%
        slice_min(
          order_by = .data[[x_col]],
          n = top_processing_proteins,
          with_ties = FALSE
        )
    )
  }

  ggplot(data) +
    aes(
      x     = .data[[x_col]],
      y     = beta,
      size  = -log10(p.value.adj),
      color = -log10(p.value.adj)
    ) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "red", linewidth = 0.4) +
    scale_x_continuous(limits = ggpmisc::symmetric_limits) +
    scale_y_continuous(limits = ggpmisc::symmetric_limits) +
    scale_color_viridis_c(option = "plasma") +
    theme_bw() +
    labs(
      title  = title_label,
      x      = paste0("Processing log\u2082 FC\n(", x_col, ")"),
      y      = "Cluster \u03b2 (Small vs. Big)",
      color  = "-log\u2081\u2080 adj. P",
      size   = "-log\u2081\u2080 adj. P"
    )
}

# Map each condition to the appropriate effect-size column.
# Column names are taken from the SomaLogic CSVs; adjust if yours differ.
condition_map <- list(
  FedFasted    = list(data = sens_prot_list$FedFasted,    x = "Time_12_hours", label = "Fed-Fasted State"),
  FreezeThaw   = list(data = sens_prot_list$FreezeThaw,   x = "cycle_10",      label = "Freeze-Thaw Cycle"),
  TimeToDecant = list(data = sens_prot_list$TimeToDecant, x = "Time_24_hours", label = "Time to Decant"),
  TimeToFreeze = list(data = sens_prot_list$TimeToFreeze, x = "Time_24_hours", label = "Time to Freeze"),
  TimeToSpin   = list(data = sens_prot_list$TimeToSpin,   x = "Time_24_hours", label = "Time to Spin")
)

# Keep only conditions for which we actually have data
condition_map <- Filter(function(cond) !is.null(cond$data), condition_map)

panels <- lapply(condition_map, function(cond) {
  make_diag_panel(cond$data, cond$x, cond$label, TOP_PROCESSING_PROTEINS)
})

diag_figure <- ggpubr::ggarrange(
  plotlist      = panels,
  common.legend = TRUE,
  legend        = "top"
)
# diag_figure is returned to the calling environment.
# index.qmd calls print(diag_figure) to render it at the correct position.
