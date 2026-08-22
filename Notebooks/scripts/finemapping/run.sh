#!/usr/bin/env bash
# Purpose: Render the fine-mapping exercise notebook to $HOME/fm_out
# Author:  course materials
# Date:    2026-06-03
# Usage:   bash run.sh   (run from session-finemapping/)
set -euo pipefail

out_dir="${HOME}/fm_out"
mkdir -p "${out_dir}"

echo "Rendering 01_finemapping.qmd → ${out_dir}/"
quarto render 01_finemapping.qmd --output-dir "${out_dir}"
echo "Done: ${out_dir}/01_finemapping.html"
