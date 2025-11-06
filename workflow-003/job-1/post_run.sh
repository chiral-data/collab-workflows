#!/bin/bash
set -e

mkdir outputs/
cp -r ./ligand_library ./outputs/
cp ./config.txt ./outputs/
cp -r ./4OHU_A_NAD_fixed_with_NAD.pdb ./outputs/
cp -r ./results ./outputs/
cp ./variants.svg ./outputs/
