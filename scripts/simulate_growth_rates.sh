#!/bin/bash

python3 simulate_growth_rates.py \
    --subject_id F01 \
    --model_name agora103_genus.qza \
    --pickled_gsmm_out ../data/pickled_models/pickled_f01_gurobi_2 \
    --solver gurobi \
    --threads 10 \
    --diet_fp ../data/diets/vmh_eu_average_agora.qza \
    --tradeoff 0.7 \
    --growth_out_fp ../data/growth_rates/growth_f01_gurobi_euavg_07.zip 