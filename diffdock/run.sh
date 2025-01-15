#!/bin/bash
# python3.8 -m inference --config default_inference_args.yaml  --protein_ligand_csv data/splittt/siema1.csv --out_dir results/ipz5
# python3.8 -m inference --config default_inference_args.yaml  --protein_ligand_csv data/splittt/siema2.csv --out_dir results/ipz6

# python3.8 -m inference --config default_inference_args.yaml  --protein_ligand_csv data/nowy_diffdock/nowy_dock0.csv --out_dir results/nowy_diffdock0
# python3.8 -m inference --config default_inference_args.yaml  --protein_ligand_csv data/nowy_diffdock/nowy_dock1000.csv --out_dir results/nowy_diffdock1000
# python3.8 -m inference --config default_inference_args.yaml  --protein_ligand_csv data/nowy_diffdock/nowy_dock2000.csv --out_dir results/nowy_diffdock2000

python3.8 -m inference --config default_inference_args.yaml  --protein_ligand_csv data/inhibitors.csv --out_dir results/inhibitors