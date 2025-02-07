import micom
import os
import pandas as pd
import shutil
import logging

# Enable verbose MICOM logging
micom.logger.setLevel(logging.DEBUG)

# Paths and constants
QZA_DIR = "../data/qiime_outputs/"
SUBJECT_IDS = ["F01", "M01", "M02"]
MODEL_DIR = "../data/models/"
MODEL_FP = os.path.join(MODEL_DIR, "agora103_genus.qza")
MODEL_EXTRACT_FP = os.path.join(MODEL_DIR, "agora103_genus")
PICKLED_MODELS = "../data/pickled_models/" #add separate subdirectories for each subject?
GROWTH_OUTPUT = "../data/growth_rates/" #add separate subdirectories for each subject here, too?

# Ensure clean directories (do I care about this - I don't think I need to re-run models once we figure this out)
#if os.path.exists(PICKLED_MODELS):
    #shutil.rmtree(PICKLED_MODELS)  # Clear previous models
os.makedirs(PICKLED_MODELS, exist_ok=True)
os.makedirs(GROWTH_OUTPUT, exist_ok=True)

def load_subject_data(subject_id, qza_dir, collapse_on='genus'):
    """
    Load subject data and convert it to a MICOM-compatible format.
    """
    feature_table_fp = os.path.join(qza_dir, f"{subject_id}_feature_table.qza")
    taxonomy_fp = os.path.join(qza_dir, f"{subject_id}_taxonomy.qza")
    subject_micom = micom.taxonomy.qiime_to_micom(
        feature_table_fp, taxonomy_fp, collapse_on=collapse_on)
    return subject_micom

# Load data for F01 and filter for a single timepoint
f01_micom = load_subject_data("F01", QZA_DIR)
#script ran properly on single timepoint (f01_micom_test)
#f01_micom_test = f01_micom[f01_micom["sample_id"] == "1201737600"]

# Load AGORA genus database
agora_genus_db = micom.qiime_formats.load_qiime_model_db(MODEL_FP, MODEL_EXTRACT_FP)

# Build MICOM community models
# took 2.5h to run on all F01 timepoints with 1 thread
from micom.workflows import build
manifest = build(
    f01_micom,
    out_folder=PICKLED_MODELS, # can I name these f"{subject_id}_" + sample_id(epoch_time)?
    model_db=MODEL_FP,
    solver="gurobi",
    threads=1  #Getting an error when I try to run this with >1 thread
)
print(manifest)

# Load the medium for simulations
from micom.qiime_formats import load_qiime_medium
western_diet = load_qiime_medium("../data/diets/western_diet_gut_agora.qza")

# Simulate growth rates
# printing progress bar and reading time + row, column, nonzero values for each model separately
from micom.workflows import grow, save_results
growth = grow(
    manifest,
    PICKLED_MODELS,
    medium=western_diet,
    tradeoff=0.7,
    threads=5  
)

# Save results
save_results(growth, os.path.join(GROWTH_OUTPUT, "growth.zip"))
print(f"Growth simulation completed. Results saved to {GROWTH_OUTPUT}")