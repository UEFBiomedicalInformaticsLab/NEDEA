#!/usr/bin/env bash


if [[ ! -d "OutPutFiles/Compare_external/drug2ways/" ]]; then
  mkdir -p "OutPutFiles/Compare_external/drug2ways/"
fi

python -m drug2ways explore \
       --graph=InputFiles/Compare_external/drug2ways/drug2ways_network.tsv \
       --fmt=tsv \
       --sources=InputFiles/Compare_external/drug2ways/drug2ways_drugs.tsv \
       --targets=InputFiles/Compare_external/drug2ways/drug2ways_phenotypes.tsv \
       --output=OutPutFiles/Compare_external/drug2ways \
       --name test_job \
       --lmax=4