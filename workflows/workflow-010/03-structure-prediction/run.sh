#!/bin/bash
set -e
python3 node3.py --inputs inputs --outputs outputs --checkpoint /workflow/abodybuilder3/output/plddt-loss/best_second_stage.ckpt --device "${PARAM_DEVICE:-cuda}"
