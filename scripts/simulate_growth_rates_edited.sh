#!/bin/bash

python3 simulate_growth_rates_edited.py \
    --subject_id F01 \
    --model_name agora103_genus.qza \
    --pickled_gsmm_out ../data/pickled_models/pickled_f01_gurobi \
    --solver gurobi \
    --threads 10 \
    --diet_fp ../data/diets/vmh_high_fiber_agora.qza \
    --tradeoff 0.3 \
    --growth_out_fp ../data/growth_rates/growth_f01_gurobi_vmhfiber_03.zip 