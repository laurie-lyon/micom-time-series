#!/bin/bash

python3 simulate_growth_rates_edited.py \
    --subject_id M02 \
    --model_name agora201_refseq216_genus_1.qza \
    --pickled_gsmm_out ../data/pickled_models/pickled_M02_agora201_gurobi \
    --solver gurobi \
    --threads 10 \
    --diet_fp ../data/diets/western_diet_gut_agora.qza \
    --tradeoff 0.3 \
    --growth_out_fp ../data/growth_rates/growth_M02_agora201_gurobi_wd_03.zip 
