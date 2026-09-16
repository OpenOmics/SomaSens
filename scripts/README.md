# Scripts

This directory contains the source code for the SomaSens quality-control and
processing-sensitive protein scoring workflow.

## Main Workflow

- `00_prepare_data.r`: read the local SomaScan ADAT file, merge phenotype
  metadata, harmonize fields used downstream, and save `data/adata_meta.rds`.
- `01_cluster_diagnosis.r`: compare unexplained sample clusters against
  SomaLogic processing-sensitive protein signatures.
- `02_fedfast_scoring.r`: calculate fed-fasted scores using both log2 fold
  change thresholding and independent protein selection.
- `03_de_analysis.r`: run differential-expression models after score-based
  filtering, including the Dream split-eBayes workflow used in the report.
- `index.qmd`: Quarto source for the detailed pre-analytical variation page and
  the standalone collaborator report. The website landing page and other topic
  pages live at the repository root.

The older development report and exploratory scripts are not required for the
focused GitHub workflow. The reusable pieces have been moved into the scripts
listed above.
