# =============================================================================
# 00_prepare_data.R
#
# Purpose: Load and prepare the SomaScan dataset for downstream analysis.
#          Produces `adata_meta`, a log2-transformed ADAT joined with sample
#          metadata and augmented with a k-means cluster assignment column
#          ('Big' / 'Small').
#
# This script is sourced automatically by 01_cluster_diagnosis.r when no
# cached de_cluster_results.rds exists. It can also be run interactively.
#
# Inputs (all paths relative to the repo root; adjust Section 0 as needed):
#   data/
#     <adat_file>        - fully normalised SomaScan ADAT
#     <phenotype_file>   - sample phenotype / metadata (Excel or CSV)
#     <shipment_file>    - optional shipment/wave information (Excel), or NULL
#
# Output:
#   adata_meta           - object placed in the calling environment; also
#                          cached to data/adata_meta.rds for faster re-runs
#
# Dependencies: tidyverse, SomaDataIO, readxl, here
# =============================================================================

library(tidyverse)
library(SomaDataIO)
library(here)

# =============================================================================
# 0. User-configurable paths and settings
# =============================================================================

ADAT_FILE     <- here("data", "CHI-24-013_v5.0_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
PHENO_FILE    <- here("data", "Final_486_phenotype_merged_07242025.xlsx")
SHIPMENT_FILE <- here("data", "CHI-24-013_Shipment_Wave_Information_26AUG2025_TL.xlsx")  # NULL if not available

# Rename map applied to phenotype table: c(new_name = "old_name")
PHENO_RENAME <- c(SampleId = "SID", Sex = "Reported Sex")

# SampleId patterns to exclude after joining
EXCLUDE_SAMPLE_PATTERNS <- "SML_QC2_220143"

# Duplicate samples to remove, identified by (SampleId, PlatePosition)
DUPLICATE_REMOVALS <- list(
  list("10202019", "C8"),
  list("10207007", "D12"),
  list("10205093", "H10"),
  list("10207090", "B11")
)

# k-means settings
KMEANS_K    <- 2
KMEANS_SEED <- 100
KMEANS_ALG  <- "MacQueen"

# Cache path — avoids re-running everything on subsequent renders
ADATA_CACHE <- here("data", "adata_meta.rds")


# =============================================================================
# 1. Return early if cache exists
# =============================================================================

if (file.exists(ADATA_CACHE)) {
  message("Loading cached adata_meta from ", ADATA_CACHE)
  adata_meta <<- readRDS(ADATA_CACHE)
  return(invisible(NULL))
}


# =============================================================================
# 2. Read ADAT
# =============================================================================

if (!file.exists(ADAT_FILE)) {
  stop(
    "ADAT file not found:\n  ", ADAT_FILE, "\n",
    "Place your normalised .adat file in data/ and update ADAT_FILE in\n",
    "scripts/00_prepare_data.r."
  )
}

message("Reading ADAT ...")
adata <- SomaDataIO::read_adat(ADAT_FILE)
adata <- adata %>%
  mutate(SampleType = ifelse(grepl("QC", SampleId), "QC2", SampleType))


# =============================================================================
# 3. Read and join phenotype metadata
# =============================================================================

message("Reading phenotype data ...")
pheno <- if (grepl("\\.xlsx?$", PHENO_FILE, ignore.case = TRUE)) {
  readxl::read_excel(PHENO_FILE)
} else {
  read.csv(PHENO_FILE)
}

pheno <- pheno %>%
  rename(any_of(PHENO_RENAME)) %>%
  mutate(
    SampleId = as.character(SampleId),
    across(any_of(c("full_PC1", "full_PC2")), as.numeric)
  )

if (!is.null(SHIPMENT_FILE) && file.exists(SHIPMENT_FILE)) {
  message("Reading shipment data ...")
  ship <- readxl::read_excel(SHIPMENT_FILE) %>%
    rename(any_of(c(SampleId = "Final SID", SampleId = "SampleId"))) %>%
    select(-any_of("Barcode")) %>%
    distinct()
  pheno <- left_join(pheno, ship, by = "SampleId")
}

adata_meta <- left_join(adata, pheno, by = "SampleId")


# =============================================================================
# 4. Remove known duplicate / low-quality samples
# =============================================================================

for (dup in DUPLICATE_REMOVALS) {
  adata_meta <- adata_meta[
    !(adata_meta$SampleId == dup[[1]] & adata_meta$PlatePosition == dup[[2]]), ]
}


# =============================================================================
# 5. Filter to biological samples and log2-transform
# =============================================================================

message("Filtering and log2-transforming ...")
adata_meta <- adata_meta %>%
  filter(
    SampleType == "Sample",
    !str_detect(SampleId, EXCLUDE_SAMPLE_PATTERNS)
  ) %>%
  log2()

rownames(adata_meta) <- adata_meta$SampleId


# =============================================================================
# 6. PCA + k-means clustering
# =============================================================================
# Cluster on scaled log2 RFU of human analytes.
# The smaller cluster is labelled "Small"; the larger "Big".

message("Running k-means clustering ...")

human_apts <- getAnalyteInfo(adata_meta) %>%
  filter(Organism == "Human") %>%
  pull(AptName)

adata_scaled <- scale(as.matrix(adata_meta[, human_apts]))

set.seed(KMEANS_SEED)
km <- kmeans(adata_scaled, centers = KMEANS_K, algorithm = KMEANS_ALG)

small_label       <- as.integer(names(which.min(table(km$cluster))))
adata_meta$Cluster <- ifelse(km$cluster == small_label, "Small", "Big")

message(sprintf(
  "Cluster sizes:  Big = %d   Small = %d",
  sum(adata_meta$Cluster == "Big"),
  sum(adata_meta$Cluster == "Small")
))


# =============================================================================
# 7. Cache and export to calling environment
# =============================================================================

saveRDS(adata_meta, ADATA_CACHE)
message("adata_meta saved to ", ADATA_CACHE)

# Make adata_meta available in the environment that sourced this script
adata_meta <<- adata_meta
