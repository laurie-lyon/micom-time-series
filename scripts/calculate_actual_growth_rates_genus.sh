#!/bin/bash

python3 ./calculate_actual_growth_rates_genus.py \
    --feature_tables_dir ../data/qiime_outputs/ \
    --taxonomy_dir ../data/qiime_outputs \
    --output_dir ../data/actual_growth_rates_genus_percent/