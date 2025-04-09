#!/usr/bin/env bash

## This script is being used to test the new functionalities of simulate_growth_rates_edited.py
## e.g. the added_metabolites and manifest_summary output files
## I also want to see if redoing the simulations with no other changes produces different results
## e.g. will the same parameter combos be ranked the same way and will the same genera be best predicted 

# Load diet shorthand mappings from config file
source diet_config.sh

# Define arrays of parameters
SUBJECT_IDS=("M01")
DIETS=("vmh_eu_average_agora.qza" "western_diet_gut_agora.qza" "vmh_high_fiber_agora.qza" "vmh_high_fat_low_carb_agora.qza")

# Define tradeoff values explicitly to prevent "00" issue
TRADEOFFS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

# Set paths
MODEL_NAME="agora201_refseq216_genus_1.qza"
PICKLED_DIR="../data/pickled_models/"
DIET_DIR="../data/diets/"
GROWTH_OUT_DIR="../data/growth_rates/"
SOLVER="gurobi"
THREADS=10

# Function to get diet shorthand using environment variables
get_diet_shorthand() {
    local diet_key="DIET_SHORTHAND_${1%.qza}"
    local shorthand="${!diet_key}"
    echo "${shorthand:-${1%.qza}}"  # Fallback: Remove .qza if not found
}

# Loop through all combinations
for SUBJECT_ID in "${SUBJECT_IDS[@]}"; do
    for DIET in "${DIETS[@]}"; do
        for TRADEOFF in "${TRADEOFFS[@]}"; do
            
            # Convert tradeoff to shorthand (0.1 → 01, 1.0 → 10)
            TRADEOFF_SHORT=$(printf "%02d" "$(echo "$TRADEOFF * 10 / 1" | bc)")

            # Get diet shorthand
            DIET_SHORT=$(get_diet_shorthand "$DIET")

            # Generate output filenames dynamically
            #added a 2 at the end of these folder names to specify a 2nd run with same parameters as before 
            PICKLED_OUT="${PICKLED_DIR}pickled_${SUBJECT_ID}_agora201_gurobi_2"
            GROWTH_OUT="${GROWTH_OUT_DIR}growth_${SUBJECT_ID}_agora201_gurobi_${DIET_SHORT}_${TRADEOFF_SHORT}_2.zip"

            # Run the simulation
            echo "Running simulation for Subject: $SUBJECT_ID, Diet: $DIET_SHORT, Tradeoff: $TRADEOFF_SHORT"

            python3 simulate_growth_rates_edited.py \
                --subject_id "$SUBJECT_ID" \
                --model_name "$MODEL_NAME" \
                --pickled_gsmm_out "$PICKLED_OUT" \
                --solver "$SOLVER" \
                --threads "$THREADS" \
                --diet_fp "${DIET_DIR}${DIET}" \
                --tradeoff "$TRADEOFF" \
                --growth_out_fp "$GROWTH_OUT"

            echo "Completed: $GROWTH_OUT"
            echo "--------------------------------------"
        done
    done
done