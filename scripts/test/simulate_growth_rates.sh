#!/bin/bash

python3 simulate_growth_rates.py \
    --subject_id F01 \
    --model_name agora103_genus.qza \
    --pickled_gsmm_out ../../data/pickled_models/test \
    --solver gurobi \
    --threads 10 \
    --diet_fp ../../data/diets/western_diet_gut_agora.qza \
    --tradeoff 0.3 \
    --growth_out_fp ../../data/growth_rates/test_growth.zip 