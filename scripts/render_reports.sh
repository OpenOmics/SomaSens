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

if [[ "${target}" == "all" || "${target}" == "docs" ]]; then
  echo "Rendering GitHub Pages site: ${repo_dir}/docs/"
  cd "${repo_dir}"
  quarto render --cache-refresh
fi

if [[ "${target}" == "all" || "${target}" == "standalone" ]]; then
  if [[ "${target}" == "standalone" ]]; then
    echo "Rendering detailed report for standalone export"
    cd "${repo_dir}"
    quarto render scripts/index.qmd --cache-refresh
  fi

  echo "Creating standalone collaborator report: ${repo_dir}/SomaSens_standalone.html"
  cp "${repo_dir}/docs/scripts/index.html" "${repo_dir}/SomaSens_standalone.html"
fi

echo "Done."
