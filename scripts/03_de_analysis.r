# =============================================================================
# 03_de_analysis.R
#
# Purpose: Reusable differential-expression helpers for SomaScan RFU data.
#
# This script is a cleaned, repo-focused version of the model helpers in
# scripts/old_scripts/. It keeps the same core ideas:
#   - limma for fixed-effect linear modeling on log2 SomaScan RFU
#   - dream for random-effect modeling with site as a random intercept
#   - optional lmerTest/lmerSeq entry points for model comparison
#   - split eBayes by SomaScan dilution group when requested
#   - volcano plotting for model output
#
# The example runs used by index.qmd fit the same Dream split-eBayes model
# before and after filtering the higher-score GMM group from the independent
# fed-fasted score at cor < 0.75.
#
# Dependencies: tidyverse, SomaDataIO, Biobase, limma, variancePartition,
#               lmerTest, ggh4x, here
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SomaDataIO))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(here))

source(here("scripts", "02_fedfast_scoring.r"))

# =============================================================================
# 0. User-configurable settings
# =============================================================================

DE_CACHE_ALL_SAMPLES <- here("data", "de_dream_all_samples.rds")
DE_CACHE_INDEP075 <- here("data", "de_dream_indep075_filter.rds")
PHENOTYPE_COL <- "Active_Asthma"
DEFAULT_COVARIATES <- c("Age", "Sex", "full_PC1", "full_PC2")
DEFAULT_RANDOM_EFFECT <- "site"


# =============================================================================
# 1. Metadata preparation
# =============================================================================

load_prepared_adata <- function(cache_path = here("data", "adata_meta.rds")) {
  if (exists("adata_meta")) {
    return(adata_meta)
  }

  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  source(here("scripts", "00_prepare_data.r"), local = FALSE)

  if (!exists("adata_meta")) {
    stop("adata_meta could not be loaded or built.")
  }

  adata_meta
}

recode_active_asthma <- function(x) {
  case_when(
    x == "Active Asthma" ~ "Case",
    x == "Case" ~ "Case",
    x == "Control" ~ "Control",
    TRUE ~ NA_character_
  )
}

impute_pc_by_site <- function(data, pc_cols = c("full_PC1", "full_PC2"), site_col = "site") {
  for (pc_col in intersect(pc_cols, colnames(data))) {
    if (!any(is.na(data[[pc_col]]))) {
      next
    }

    if (site_col %in% colnames(data)) {
      data <- data %>%
        group_by(.data[[site_col]]) %>%
        mutate(
          "{pc_col}" := ifelse(
            is.na(.data[[pc_col]]),
            mean(.data[[pc_col]], na.rm = TRUE),
            .data[[pc_col]]
          )
        ) %>%
        ungroup()
    }

    if (any(is.na(data[[pc_col]]))) {
      data[[pc_col]][is.na(data[[pc_col]])] <- mean(data[[pc_col]], na.rm = TRUE)
    }
  }

  data
}

prepare_model_input <- function(adata_meta,
                                phenotype_col = PHENOTYPE_COL,
                                covariates = DEFAULT_COVARIATES,
                                random_effect = NULL,
                                human_only = TRUE) {
  model_data <- adata_meta

  if (phenotype_col == "Active_Asthma") {
    model_data[[phenotype_col]] <- recode_active_asthma(model_data[[phenotype_col]])
  }

  model_data <- impute_pc_by_site(model_data)

  needed_cols <- unique(c("SampleId", phenotype_col, covariates, random_effect))
  missing_cols <- setdiff(needed_cols, colnames(model_data))
  if (length(missing_cols) > 0) {
    stop("Missing model metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  model_data <- model_data %>%
    filter(if_all(all_of(setdiff(needed_cols, "SampleId")), ~ !is.na(.x)))

  if (is.character(model_data[[phenotype_col]])) {
    phenotype_values <- sort(unique(model_data[[phenotype_col]]))
    if (all(c("Control", "Case") %in% phenotype_values)) {
      model_data[[phenotype_col]] <- factor(model_data[[phenotype_col]], levels = c("Control", "Case"))
    } else {
      model_data[[phenotype_col]] <- factor(model_data[[phenotype_col]])
    }
  }

  for (col in intersect(c("Sex", random_effect), colnames(model_data))) {
    model_data[[col]] <- factor(model_data[[col]])
  }

  adata_es <- SomaDataIO::adat2eSet(model_data)
  analytes <- SomaDataIO::getAnalyteInfo(model_data) %>%
    select(AptName, SeqId, Target, UniProt, EntrezGeneID, Dilution, Organism)

  if (human_only) {
    analytes <- analytes %>%
      filter(Organism == "Human") %>%
      filter(Dilution %in% c(20, 0.5, 0.005))
  }

  analytes <- analytes %>% filter(AptName %in% rownames(Biobase::exprs(adata_es)))
  expr_mat <- Biobase::exprs(adata_es)[analytes$AptName, , drop = FALSE]
  meta <- Biobase::pData(adata_es)

  list(
    data = model_data,
    expr = expr_mat,
    meta = meta,
    analytes = analytes
  )
}


# =============================================================================
# 2. Shared model output helpers
# =============================================================================

get_contrast_name <- function(phenotype_col, phenotype_values) {
  if (is.numeric(phenotype_values)) {
    return(phenotype_col)
  }

  if (all(c("Control", "Case") %in% levels(phenotype_values))) {
    return(paste0(phenotype_col, "Case-", phenotype_col, "Control"))
  }

  stop("Only numeric phenotypes or Control/Case factors are currently supported.")
}

top_table_by_dilution <- function(fit,
                                  analytes,
                                  coef,
                                  split_eBayes = TRUE,
                                  p_adjust_method = "BH") {
  if (!split_eBayes) {
    fit_ebayes <- limma::eBayes(fit)
    return(
      limma::topTable(
        fit_ebayes,
        coef = coef,
        number = Inf,
        sort.by = "none",
        adjust.method = p_adjust_method
      ) %>%
        tibble::rownames_to_column("AptName") %>%
        left_join(analytes, by = "AptName")
    )
  }

  split(analytes, analytes$Dilution) %>%
    map_dfr(function(analyte_group) {
      rows <- intersect(analyte_group$AptName, rownames(fit))
      if (length(rows) == 0) {
        return(tibble())
      }

      fit_ebayes <- limma::eBayes(fit[rows, ])
      limma::topTable(
        fit_ebayes,
        coef = coef,
        number = Inf,
        sort.by = "none",
        adjust.method = p_adjust_method
      ) %>%
        tibble::rownames_to_column("AptName") %>%
        left_join(analytes, by = "AptName")
    })
}


# =============================================================================
# 3. Model runners
# =============================================================================

fit_limma_model <- function(adata_meta,
                            phenotype_col = PHENOTYPE_COL,
                            covariates = c(DEFAULT_RANDOM_EFFECT, DEFAULT_COVARIATES),
                            split_eBayes = TRUE) {
  model_input <- prepare_model_input(
    adata_meta = adata_meta,
    phenotype_col = phenotype_col,
    covariates = covariates,
    random_effect = NULL
  )

  covariate_str <- paste(covariates, collapse = " + ")
  design_formula <- as.formula(paste0("~ 0 + ", phenotype_col, " + ", covariate_str))
  design <- stats::model.matrix(design_formula, data = model_input$meta)

  if (is.numeric(model_input$meta[[phenotype_col]])) {
    fit <- limma::lmFit(model_input$expr, design = design)
    coef_name <- phenotype_col
  } else {
    contrast_name <- get_contrast_name(phenotype_col, model_input$meta[[phenotype_col]])
    contrasts <- limma::makeContrasts(contrasts = contrast_name, levels = design)
    fit <- limma::lmFit(model_input$expr, design = design) %>%
      limma::contrasts.fit(contrasts)
    coef_name <- contrast_name
  }

  top_table_by_dilution(
    fit = fit,
    analytes = model_input$analytes,
    coef = coef_name,
    split_eBayes = split_eBayes
  )
}

fit_dream_model <- function(adata_meta,
                            phenotype_col = PHENOTYPE_COL,
                            covariates = DEFAULT_COVARIATES,
                            random_effect = DEFAULT_RANDOM_EFFECT,
                            covar_score = NULL,
                            split_eBayes = TRUE) {
  model_input <- prepare_model_input(
    adata_meta = adata_meta,
    phenotype_col = phenotype_col,
    covariates = c(covariates, covar_score),
    random_effect = random_effect
  )

  fixed_terms <- paste(c(phenotype_col, covariates, covar_score), collapse = " + ")
  formula <- as.formula(paste0("~ 0 + ", fixed_terms, " + (1|", random_effect, ")"))

  if (is.numeric(model_input$meta[[phenotype_col]])) {
    dream_fit <- variancePartition::dream(
      exprObj = model_input$expr,
      formula = formula,
      data = model_input$meta
    )
    coef_name <- phenotype_col
  } else {
    contrast_name <- paste0(phenotype_col, "Case - ", phenotype_col, "Control")
    contrasts <- variancePartition::makeContrastsDream(
      formula,
      model_input$meta,
      contrasts = c(Case_vs_Control = contrast_name)
    )
    dream_fit <- variancePartition::dream(
      exprObj = model_input$expr,
      formula = formula,
      data = model_input$meta,
      L = contrasts
    )
    coef_name <- "Case_vs_Control"
  }

  top_table_by_dilution(
    fit = dream_fit,
    analytes = model_input$analytes,
    coef = coef_name,
    split_eBayes = split_eBayes
  )
}

fit_lmer_test_model <- function(adata_meta,
                                phenotype_col = PHENOTYPE_COL,
                                covariates = DEFAULT_COVARIATES,
                                random_effect = DEFAULT_RANDOM_EFFECT) {
  model_input <- prepare_model_input(
    adata_meta = adata_meta,
    phenotype_col = phenotype_col,
    covariates = covariates,
    random_effect = random_effect
  )

  if (!all(c("Control", "Case") %in% levels(model_input$meta[[phenotype_col]]))) {
    stop("fit_lmer_test_model currently expects a Control/Case phenotype.")
  }

  meta <- model_input$meta
  fixed_terms <- paste(c(phenotype_col, covariates), collapse = " + ")

  map_dfr(model_input$analytes$AptName, function(aptn) {
    df <- meta
    df$Protein <- as.numeric(model_input$expr[aptn, ])
    form <- as.formula(paste0("Protein ~ ", fixed_terms, " + (1|", random_effect, ")"))
    model <- lmerTest::lmer(form, data = df, REML = FALSE)
    coef_tab <- summary(model)$coefficients %>% as.data.frame()
    coef_name <- paste0(phenotype_col, "Case")

    tibble(
      AptName = aptn,
      Estimate = coef_tab[coef_name, "Estimate"],
      se = coef_tab[coef_name, "Std. Error"],
      t = coef_tab[coef_name, "t value"],
      p.value = coef_tab[coef_name, "Pr(>|t|)"]
    )
  }) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    left_join(model_input$analytes, by = "AptName")
}

fit_lmer_seq_model <- function(...) {
  if (!requireNamespace("lmerSeq", quietly = TRUE)) {
    stop(
      "lmerSeq is not installed in this environment. ",
      "The original report tested it, but this focused repo treats it as optional."
    )
  }

  stop(
    "lmerSeq support is intentionally left as an optional extension. ",
    "Use fit_lmer_test_model() or fit_dream_model() for the reproducible workflow here."
  )
}


# =============================================================================
# 4. Fed-fasted independent-score filtering for the worked example
# =============================================================================

add_indep075_score_filter <- function(adata_meta, cor_threshold = 0.75) {
  score_df <- score_from_independent_proteins(
    cor_threshold = cor_threshold,
    data = adata_meta,
    score_name = "score_indep075"
  )

  gmm_fit <- fit_score_gmm(score_df, score_col = "score_indep075")

  score_data <- gmm_fit$gmm_data %>%
    select(SampleId, score_indep075, GMMClass) %>%
    rename(gmm_class_indep075 = GMMClass)

  output <- adata_meta
  match_idx <- match(output$SampleId, score_data$SampleId)
  output$score_indep075 <- score_data$score_indep075[match_idx]
  output$gmm_class_indep075 <- score_data$gmm_class_indep075[match_idx]
  output
}

filter_indep075_high_score_samples <- function(adata_meta, cor_threshold = 0.75) {
  add_indep075_score_filter(adata_meta, cor_threshold = cor_threshold) %>%
    filter(gmm_class_indep075 != "Higher score")
}

run_dream_all_samples <- function(cache_path = DE_CACHE_ALL_SAMPLES,
                                  force = FALSE) {
  if (file.exists(cache_path) && !force) {
    return(readRDS(cache_path))
  }

  prepared_adata <- load_prepared_adata()

  dream_results <- fit_dream_model(
    adata_meta = prepared_adata,
    phenotype_col = PHENOTYPE_COL,
    covariates = DEFAULT_COVARIATES,
    random_effect = DEFAULT_RANDOM_EFFECT,
    split_eBayes = TRUE
  )

  output <- list(
    method = "Dream split eBayes",
    filter = "No fed-fasted filtering; all samples retained",
    results = dream_results,
    n_samples = nrow(prepared_adata)
  )

  saveRDS(output, cache_path)
  output
}

run_dream_indep075_filter <- function(cache_path = DE_CACHE_INDEP075,
                                      force = FALSE) {
  if (file.exists(cache_path) && !force) {
    return(readRDS(cache_path))
  }

  prepared_adata <- load_prepared_adata()
  filtered_adata <- filter_indep075_high_score_samples(prepared_adata, cor_threshold = 0.75)

  dream_results <- fit_dream_model(
    adata_meta = filtered_adata,
    phenotype_col = PHENOTYPE_COL,
    covariates = DEFAULT_COVARIATES,
    random_effect = DEFAULT_RANDOM_EFFECT,
    split_eBayes = TRUE
  )

  output <- list(
    method = "Dream split eBayes",
    filter = "Independent fed-fasted score, cor < 0.75; higher-score GMM group removed",
    results = dream_results,
    n_samples_before = nrow(prepared_adata),
    n_samples_after = nrow(filtered_adata),
    gmm_counts = add_indep075_score_filter(prepared_adata, 0.75) %>%
      count(gmm_class_indep075, Cluster, site, name = "n")
  )

  saveRDS(output, cache_path)
  output
}


# =============================================================================
# 5. Volcano plot
# =============================================================================

volcano_break_axis <- function(data,
                               x = "logFC",
                               y = "adj.P.Val",
                               y_breaks = c(3, 17),
                               x_breaks = c(0.75, 1.75),
                               pCutoff = 0.1,
                               FCcutoff = 0.25,
                               title = NULL,
                               subtitle = NULL,
                               caption = NULL,
                               add.hline = -log10(0.05),
                               add.vline = NULL) {
  plot_data <- data %>%
    filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>%
    mutate(
      p_for_plot = pmax(.data[[y]], .Machine$double.xmin),
      Significance = case_when(
        abs(.data[[x]]) >= FCcutoff & .data[[y]] < pCutoff ~ "p-value and log_2FC",
        abs(.data[[x]]) >= FCcutoff ~ "Log2 FC",
        .data[[y]] < pCutoff ~ "p-value",
        TRUE ~ "NS"
      ),
      Intercept = -log10(pCutoff),
      XIntercept = ifelse(.data[[x]] < 0, -FCcutoff, FCcutoff),
      XBreak = "x",
      YBreak = "y"
    )

  if (!is.null(y_breaks)) {
    plot_data <- plot_data %>%
      mutate(
        YBreak = cut(-log10(p_for_plot), c(-Inf, y_breaks, Inf)) %>%
          forcats::fct_rev(),
        Intercept = ifelse(
          YBreak == levels(YBreak)[length(levels(YBreak))],
          -log10(pCutoff),
          NA_real_
        )
      )
  }

  if (!is.null(x_breaks)) {
    plot_data <- plot_data %>%
      mutate(
        XBreak = cut(.data[[x]], c(-Inf, x_breaks, Inf)),
        XIntercept = ifelse(
          XBreak == levels(XBreak)[1],
          ifelse(.data[[x]] < 0, -FCcutoff, FCcutoff),
          NA_real_
        )
      )
  }

  y_panel_count <- if (!is.null(y_breaks)) {
    nlevels(droplevels(plot_data$YBreak))
  } else {
    1
  }
  y_panel_sizes <- case_when(
    y_panel_count == 1 ~ list(5),
    y_panel_count == 2 ~ list(c(1, 5)),
    y_panel_count == 3 ~ list(c(1, 1, 5)),
    TRUE ~ list(rep(1, y_panel_count))
  )[[1]]

  x_panel_count <- if (!is.null(x_breaks)) {
    nlevels(droplevels(plot_data$XBreak))
  } else {
    1
  }
  x_panel_sizes <- case_when(
    x_panel_count == 1 ~ list(5),
    x_panel_count == 2 ~ list(c(5, 1)),
    x_panel_count == 3 ~ list(c(5, 1, 1)),
    TRUE ~ list(rep(1, x_panel_count))
  )[[1]]

  ggplot(plot_data, aes(x = .data[[x]], y = -log10(p_for_plot))) +
    geom_point(aes(color = Significance)) +
    geom_hline(aes(yintercept = Intercept), linetype = "dashed", na.rm = TRUE) +
    geom_vline(aes(xintercept = XIntercept), linetype = "dashed", na.rm = TRUE) +
    facet_grid(YBreak ~ XBreak, scales = "free") +
    ggh4x::force_panelsizes(rows = y_panel_sizes, cols = x_panel_sizes) +
    theme_bw() +
    scale_y_continuous(breaks = c(seq(0, 25, 0.5), max(-log10(plot_data$p_for_plot)))) +
    scale_x_continuous(breaks = c(seq(-5, 5, 0.25), max(plot_data[[x]]))) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_blank(),
      legend.position = "top"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = ""
    )
}

plot_dream_volcano <- function(de_output,
                               title,
                               subtitle,
                               y_breaks = c(3, 17)) {
  volcano_break_axis(
    data = de_output$results,
    x = "logFC",
    y = "adj.P.Val",
    y_breaks = y_breaks,
    pCutoff = 0.1,
    FCcutoff = 0.25,
    title = title,
    subtitle = subtitle,
    caption = "Proteins are highlighted when adjusted p-value < 0.1 and |log2FC| >= 0.25."
  )
}

plot_all_samples_dream_volcano <- function(de_output = run_dream_all_samples()) {
  plot_dream_volcano(
    de_output = de_output,
    title = "Active Asthma DE before fed-fasted filtering",
    subtitle = "Dream split eBayes; all samples retained",
    y_breaks = c(4, 17)
  )
}

plot_indep075_dream_volcano <- function(de_output = run_dream_indep075_filter()) {
  plot_dream_volcano(
    de_output = de_output,
    title = "Active Asthma DE after fed-fasted filtering",
    subtitle = "Dream split eBayes; independent score cor < 0.75; higher-score GMM group removed"
  )
}
