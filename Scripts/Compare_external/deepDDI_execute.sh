#!/usr/bin/env bash


if [[ ! -d "OutPutFiles/Compare_external/deepDDI/" ]]; then
  mkdir -p "OutPutFiles/Compare_external/deepDDI/"
fi


python3 External_tools/kaistsystemsbiology-deepddi/run_DeepDDI.py \
        -i External_tools/kaistsystemsbiology-deepddi/examples/input_structures.txt \
        -o OutPutFiles/Compare_external/deepDDI/ 