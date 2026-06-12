# SomaSens

SomaSens is a focused workflow for diagnosing and handling pre-analytical
variation in SomaScan proteomics data. The current version concentrates on
quality-control diagnostics and scoring for processing-sensitive proteins,
with a worked example based on a multi-site asthma SomaScan study.

The rendered documentation is designed for GitHub Pages:

https://spaul-genetics.github.io/SomaSens/

## What This Repository Contains

- `scripts/00_prepare_data.r`: prepare the SomaScan ADAT object and phenotype metadata.
- `scripts/01_cluster_diagnosis.r`: diagnose unexplained sample clusters using processing-sensitive protein lists.
- `scripts/02_fedfast_scoring.r`: calculate fed-fasted scores using thresholding and independent protein selection.
- `scripts/03_de_analysis.r`: run the focused differential-expression workflow after score-based filtering.
- `scripts/index.qmd`: source document for the GitHub Pages report.
- `data/README.md`: expected input files and generated local cache files.
- `packages.R`: R package dependencies needed to reproduce the workflow.

Raw SomaScan ADAT files, phenotype spreadsheets, SomaLogic processing-sensitive
protein files, and generated `.rds` caches are intentionally not committed to
GitHub. Place those files under `data/` as described in `data/README.md`.

## Render The Documentation

From the repository root:

```bash
quarto render scripts/index.qmd --output-dir ../docs --cache-refresh
```

The command writes the GitHub Pages-ready HTML to `docs/index.html`. The
rendered HTML is committed so the site can be served directly from the
`main` branch using the `/docs` folder.

## Setup

Install the R dependencies with:

```bash
./packages.R
```

The analysis expects the local data files described in `data/README.md`.
Because these files can be large or controlled-access, they remain local and
are ignored by git.
