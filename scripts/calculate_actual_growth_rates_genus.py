#!/usr/bin/env python3
"""
Calculate Actual Growth Rates for Multiple Subjects
---------------------------------------------------

Purpose:
This script calculates the actual change in OTU abundance between consecutive samples
from QIIME2 feature tables (`feature_table.qza`) for multiple subjects. Each subject's 
data is processed independently, and the results are saved separately.

Workflow:
1. Export the feature table for each subject from QIIME2 (`.qza`) to `.biom` format.
2. Convert the `.biom` file to a readable `.tsv` format.
3. Load the `.tsv` into a pandas DataFrame and sort columns (samples) by epoch time.
4. Calculate the change in abundance (ΔAbundance) for each OTU between consecutive days.
5. Save the resulting growth rates as a `.csv` file for each subject.

Inputs:
- A directory of QIIME2 feature tables (`feature_table.qza`) for multiple subjects.
- Sample IDs (column headers) as epoch time values in seconds.

Outputs:
- One `actual_growth_rates_<subject_id>.csv` file per subject:
  - Rows = OTU IDs.
  - Columns = Epoch time of the first day in consecutive days.
  - Values = Change in abundance (ΔAbundance).

Dependencies:
- QIIME2 must be installed and accessible in the environment.
- BIOM CLI must be installed for file format conversions.
- Python libraries: pandas, numpy, subprocess, pathlib.

Usage:
    python calculate_actual_growth_rates_multi.py \
        --feature_tables_dir <path_to_feature_tables_dir> \
        --output_dir <output_directory>

Example:
    python calculate_actual_growth_rates_multi.py \
        --feature_tables_dir ../data/qiime_outputs/ \
        --output_dir ../data/actual_growth_rates/

Author: Laurie Lyon 
Date: 01/21/2025
"""

import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import argparse
import os

def load_taxonomy(taxonomy_qza, output_dir):
    """
    Load taxonomy data from a QIIME2 taxonomy artifact.
    """
    # Create a temporary directory for exporting taxonomy
    export_dir = Path(output_dir) / "taxonomy_exports"
    export_dir.mkdir(parents=True, exist_ok=True)

    # Export the taxonomy.qza file to a readable format
    subprocess.run([
        "qiime", "tools", "export",
        "--input-path", taxonomy_qza,
        "--output-path", str(export_dir)
    ], check=True)

    # Load the taxonomy TSV file
    taxonomy_tsv = export_dir / "taxonomy.tsv"
    taxonomy = pd.read_csv(taxonomy_tsv, sep="\t", index_col=0)

    # Extract genus-level taxonomy
    taxonomy["Genus"] = taxonomy["Taxon"].str.split(";").str[5]  # Extract genus-level taxonomy
    taxonomy["Genus"] = taxonomy["Genus"].str.strip()  # Clean whitespace

    return taxonomy

def process_subject(feature_table_qza, taxonomy_qza, output_dir, subject_id):
    """
    Process a single subject's feature table to calculate actual growth rates by genus.

    Parameters:
        - feature_table_qza (str): Path to the subject's feature table (.qza).
        - output_dir (str): Directory where the subject's output will be saved.
        - subject_id (str): Identifier for the subject (e.g., "F01").
        - taxonomy_qza (str): Path to the subject's taxonomy file (.qza).
    """
    # Load taxonomy data from load_taxonomy function
    taxonomy = load_taxonomy(taxonomy_qza, output_dir)

    # Create an export directory for storing intermediate files for this subject
    export_dir = Path(output_dir) / "qiime_exports" / subject_id
    export_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Export the feature table to .biom format using QIIME2
    print(f"Exporting feature table for subject {subject_id}...")
    subprocess.run([
        "qiime", "tools", "export",              # QIIME2 command to export data
        "--input-path", feature_table_qza,      # Path to the input feature table (.qza)
        "--output-path", str(export_dir)        # Directory where the .biom file will be exported
    ], check=True)

    # Path to the exported .biom file
    feature_table_biom = export_dir / "feature-table.biom"

    # Step 2: Convert the .biom file to a tab-delimited .tsv file using BIOM CLI
    print(f"Converting feature table for subject {subject_id} to TSV format...")
    feature_table_tsv = export_dir / f"{subject_id}_feature_table.tsv"

    subprocess.run([
        "biom", "convert",                      # BIOM CLI command to convert file formats
        "-i", str(feature_table_biom),          # Input .biom file
        "-o", str(feature_table_tsv),           # Output .tsv file
        "--to-tsv"                              # Convert to TSV format
    ], check=True)

    # Step 3: Load the TSV file into a pandas DataFrame for analysis
    print(f"Loading feature table for subject {subject_id}...")
    feature_table = pd.read_csv(feature_table_tsv, sep="\t", skiprows=1, index_col=0)
    #how many non-zero values are in the feature table, filter out if .1 or more are zero
    
    # Step 4: Sort sample columns by epoch time (numerically, not alphabetically)
    print(f"Sorting samples by epoch time for subject {subject_id}...")
    sorted_columns = sorted([int(col) for col in feature_table.columns])  # Convert column names to integers
    feature_table = feature_table[sorted(map(str, sorted_columns))]       # Reorder columns chronologically

    # Step 5: Calculate daily changes in abundance (ΔAbundance) for each OTU
    print(f"Normalizing feature table for subject {subject_id} to relative abundances...")
    
    #print statements to check that sum = rarefaction depth before normalization and 1 after normalization
    #print("Sum of each sample BEFORE normalization:\n", feature_table.sum(axis=0))
    feature_table = feature_table.div(feature_table.sum(axis=0), axis=1)  # Normalize to relative abundance

    #print("Sum of each sample AFTER normalization:\n", feature_table.sum(axis=0))
    #normalized_output_path = Path(output_dir) / f"{subject_id}_normalized_feature_table.csv"
    #feature_table.to_csv(normalized_output_path)

    print(f"Calculating actual growth rates for subject {subject_id}...")

    actual_growth_rates = {}

    # Loop through each pair of consecutive samples
    for i in range(len(sorted_columns) - 1):
        day1 = str(sorted_columns[i])      # Epoch time of the current day
        day2 = str(sorted_columns[i + 1])  # Epoch time of the next day

        # Calculate ΔAbundance (% change between days 1 and 2)
        epsilon = 1e-6  # Small value to avoid division by zero
        #actual_growth_rates[day1] = ((feature_table[day2] - feature_table[day1]) / (feature_table[day1] + epsilon)) * 100
        #trying log fold change to deal with outliers
        #actual_growth_rates[day1] = np.log2((feature_table[day2]) / (feature_table[day1]))
        actual_growth_rates[day1] = np.log2((feature_table[day2] + epsilon) / ((feature_table[day1]) + epsilon))

    # Convert the results into a DataFrame (rows = OTUs, columns = days)
    growth_rate_df = pd.DataFrame(actual_growth_rates)

    # Step 6: Map OTUs to their genera and collapse growth rates by genus
    print(f"Collapsing growth rates by genus for subject {subject_id}...")
    growth_rate_df["Genus"] = taxonomy.loc[growth_rate_df.index, "Genus"]
    growth_by_genus = growth_rate_df.groupby("Genus").sum() 

    # Step 7: Save the genus-level growth rates to a CSV file
    output_path = Path(output_dir) / f"{subject_id}_actual_growth_rates_by_genus.csv"
    growth_by_genus.to_csv(output_path)

    print(f"Actual growth rates saved for subject {subject_id} at: {output_path}")


def main(feature_tables_dir, taxonomy_dir, output_dir):
    """
    Main function to process feature tables for multiple subjects.

    Parameters:
        - feature_tables_dir (str): Directory containing QIIME2 feature tables (.qza files).
        - taxonomy_dir (str): Directory containing QIIME2 taxonomy files (.qza files).
        - output_dir (str): Directory where all outputs will be saved.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for file in os.listdir(feature_tables_dir):
        if file.endswith("_feature_table.qza"):
            subject_id = file.split("_")[0]
            feature_table_qza = os.path.join(feature_tables_dir, file)
            taxonomy_qza = os.path.join(taxonomy_dir, f"{subject_id}_taxonomy.qza")

            if os.path.exists(taxonomy_qza):
                process_subject(feature_table_qza, taxonomy_qza, output_dir, subject_id)
            else:
                print(f"Taxonomy file for {subject_id} not found. Skipping.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate genus-level actual growth rates for multiple subjects.")
    parser.add_argument("--feature_tables_dir", required=True, help="Directory containing QIIME2 feature tables (.qza files).")
    parser.add_argument("--taxonomy_dir", required=True, help="Directory containing QIIME2 taxonomy files (.qza files).")
    parser.add_argument("--output_dir", required=True, help="Directory where outputs will be saved.")
    args = parser.parse_args()

    main(args.feature_tables_dir, args.taxonomy_dir, args.output_dir)