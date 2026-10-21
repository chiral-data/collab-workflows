#!/bin/bash
set -e
python3 predict_structure.py \
  --inputs inputs \
  --outputs outputs \
  --checkpoint /workflow/abodybuilder3/output/plddt-loss/best_second_stage.ckpt \
  --checkpoint_lm /workflow/abodybuilder3/output/language-loss/best_second_stage.ckpt \
  --device "${PARAM_DEVICE:-cuda}"
