#!/usr/bin/env bash
set -euo pipefail

target="${1:-all}"

case "${target}" in
  all|docs|standalone)
    ;;
  *)
    echo "Usage: ./scripts/render_reports.sh [all|docs|standalone]" >&2
    exit 2
    ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"

cd "${script_dir}"

if [[ "${target}" == "all" || "${target}" == "docs" ]]; then
  echo "Rendering GitHub Pages report: ${repo_dir}/docs/index.html"
  quarto render index.qmd --output-dir ../docs --cache-refresh
fi

if [[ "${target}" == "all" || "${target}" == "standalone" ]]; then
  echo "Rendering standalone collaborator report: ${repo_dir}/SomaSens_standalone.html"
  quarto render index.qmd \
    --output SomaSens_standalone.html \
    --cache-refresh \
    -M html-math-method:mathml

  mv SomaSens_standalone.html "${repo_dir}/SomaSens_standalone.html"
fi

echo "Done."
