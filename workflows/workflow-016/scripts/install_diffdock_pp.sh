#!/bin/bash
set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/ketatam/DiffDock-PP.git}"
TARGET_DIR="${1:-${DIFFDOCK_PP_PATH:-$HOME/DiffDock-PP}}"

if [[ -d "$TARGET_DIR/.git" ]]; then
  echo "Updating existing DiffDock-PP clone at: $TARGET_DIR"
  git -C "$TARGET_DIR" pull --ff-only
else
  echo "Cloning DiffDock-PP into: $TARGET_DIR"
  git clone --depth=1 "$REPO_URL" "$TARGET_DIR"
fi

if ! command -v python >/dev/null 2>&1; then
  echo "python not found in PATH"
  exit 1
fi

echo "Installing DiffDock-PP into current Python environment"
python -m pip install --no-cache-dir \
  dill \
  biopandas \
  "e3nn<0.6" \
  wandb \
  tensorboard \
  tensorboardX

TORCH_BASE="$(python -c 'import torch; print(torch.__version__.split("+")[0])')"
TORCH_MM="$(python -c 'import torch; v=torch.__version__.split("+")[0].split("."); print(f"{v[0]}.{v[1]}.0")')"
CUDA_VER="$(python -c 'import torch; print(torch.version.cuda or "cpu")')"

if [[ "$CUDA_VER" == "cpu" || -z "$CUDA_VER" ]]; then
  CUDA_TAG="cpu"
else
  CUDA_TAG="cu${CUDA_VER//./}"
fi

PYG_URL="https://data.pyg.org/whl/torch-${TORCH_MM}+${CUDA_TAG}.html"
echo "Installing PyG extensions from: $PYG_URL"

python -m pip install --no-cache-dir \
  torch-scatter \
  torch-sparse \
  torch-cluster \
  torch-spline-conv \
  -f "$PYG_URL"

python -m pip install --no-cache-dir torch-geometric

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  ACTIVATE_DIR="$CONDA_PREFIX/etc/conda/activate.d"
  DEACTIVATE_DIR="$CONDA_PREFIX/etc/conda/deactivate.d"
  mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

  cat > "$ACTIVATE_DIR/diffdock_pp_path.sh" <<EOF
export DIFFDOCK_PP_PATH="$TARGET_DIR"
EOF

  cat > "$DEACTIVATE_DIR/diffdock_pp_path.sh" <<'EOF'
unset DIFFDOCK_PP_PATH
EOF

  export DIFFDOCK_PP_PATH="$TARGET_DIR"
  echo "Configured DIFFDOCK_PP_PATH for conda env at: $CONDA_PREFIX"
else
  if ! grep -q 'DIFFDOCK_PP_PATH=' "$HOME/.bashrc" 2>/dev/null; then
    echo "export DIFFDOCK_PP_PATH=\"$TARGET_DIR\"" >> "$HOME/.bashrc"
  fi
  export DIFFDOCK_PP_PATH="$TARGET_DIR"
  echo "Configured DIFFDOCK_PP_PATH in ~/.bashrc"
fi

echo "DiffDock-PP is ready: $TARGET_DIR"
echo "Detected main script: $TARGET_DIR/src/main_inf.py"
echo "Detected torch version: $TORCH_BASE"
