#!/bin/bash
# =============================================================================

CONTAINER="gwas-minimal.sif"

apptainer exec "${CONTAINER}" pixi run -m /opt/gwas/pixi.toml jupyter lab --ip=$(hostname) --port=$UID --allow-root