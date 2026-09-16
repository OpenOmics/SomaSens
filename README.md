# SomaSens

SomaSens is a multi-page resource for designing SomaScan studies and diagnosing
and handling pre-analytical or technical variation. It currently includes a
detailed processing-variation workflow, practical plate-randomization guidance,
and a scaffold for future CSF SOMAmer-filtering guidance.

The rendered documentation is designed for GitHub Pages:

https://spaul-genetics.github.io/SomaSens/

## What This Repository Contains

- `scripts/00_prepare_data.r`: prepare the SomaScan ADAT object and phenotype metadata.
- `scripts/01_cluster_diagnosis.r`: diagnose unexplained sample clusters using processing-sensitive protein lists.
- `scripts/02_fedfast_scoring.r`: calculate fed-fasted scores using thresholding and independent protein selection.
- `scripts/03_de_analysis.r`: run the focused differential-expression workflow after score-based filtering.
- `index.qmd`: overview and landing page.
- `scripts/index.qmd`: detailed pre-analytical variation report.
- `randomization.qmd`: plate-randomization guidance and design figures.
- `csf-somamer-filtering.qmd`: scaffold for the planned CSF filtering report.
- `_quarto.yml` and `styles.css`: website structure and presentation.
- `data/README.md`: expected input files and generated local cache files.
- `packages.R`: R package dependencies needed to reproduce the workflow.

Raw SomaScan ADAT files, phenotype spreadsheets, SomaLogic processing-sensitive
protein files, and generated `.rds` caches are intentionally not committed to
GitHub. Place those files under `data/` as described in `data/README.md`.

## Render The Documentation

Use the helper script from the repository root to render the website and the
standalone detailed report:

```bash
./scripts/render_reports.sh
```

This updates:

- `docs/`: the complete GitHub Pages website. These files can stay committed so
  the site can be served from the `main` branch using the `/docs` folder when
  GitHub Pages is enabled.
- `SomaSens_standalone.html`: a single-file HTML report for sharing privately
  with collaborators while the repository or GitHub Pages site is not public.
  It is copied from the self-contained rendered detail page and is intentionally
  ignored by git.

To render only one output:

```bash
./scripts/render_reports.sh docs
./scripts/render_reports.sh standalone
```

The standalone report is rendered with embedded resources and MathML so it can
be opened directly in a browser as a local HTML file.

## Setup

Install the R dependencies with:

```bash
./packages.R
```

The analysis expects the local data files described in `data/README.md`.
Because these files can be large or controlled-access, they remain local and
are ignored by git.
