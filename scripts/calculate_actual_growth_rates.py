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

def process_subject(feature_table_qza, output_dir, subject_id):
    """
    Process a single subject's feature table to calculate actual growth rates.

    Parameters:
        - feature_table_qza (str): Path to the subject's feature table (.qza).
        - output_dir (str): Directory where the subject's output will be saved.
        - subject_id (str): Identifier for the subject (e.g., "F01").
    """

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

    # Step 4: Sort sample columns by epoch time (numerically, not alphabetically)
    print(f"Sorting samples by epoch time for subject {subject_id}...")
    sorted_columns = sorted([int(col) for col in feature_table.columns])  # Convert column names to integers
    feature_table = feature_table[sorted(map(str, sorted_columns))]       # Reorder columns chronologically

    # Step 5: Calculate daily changes in abundance (ΔAbundance) for each OTU
    print(f"Calculating actual growth rates for subject {subject_id}...")
    actual_growth_rates = {}

    # Loop through each pair of consecutive samples
    for i in range(len(sorted_columns) - 1):
        day1 = str(sorted_columns[i])      # Epoch time of the current day
        day2 = str(sorted_columns[i + 1])  # Epoch time of the next day

        # Calculate ΔAbundance (difference in abundance between day1 and day2)
        actual_growth_rates[day1] = feature_table[day2] - feature_table[day1]

    # Convert the results into a DataFrame (rows = OTUs, columns = days)
    growth_rate_df = pd.DataFrame(actual_growth_rates)

    # Step 6: Save the actual growth rates to a CSV file
    output_path = Path(output_dir) / f"{subject_id}_actual_growth_rates.csv"
    growth_rate_df.to_csv(output_path)

    print(f"Actual growth rates saved for subject {subject_id} at: {output_path}")


def main(feature_tables_dir, output_dir):
    """
    Main function to process feature tables for multiple subjects.

    Parameters:
        - feature_tables_dir (str): Directory containing QIIME2 feature tables (.qza files).
        - output_dir (str): Directory where all outputs will be saved.
    """
    # Ensure the output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Loop through all files in the feature tables directory
    for file in os.listdir(feature_tables_dir):
        if file.endswith("_feature_table.qza"):  # Identify subject-specific feature tables by their naming pattern
            subject_id = file.split("_")[0]      # Extract subject ID from the filename (e.g., "F01")
            feature_table_qza = os.path.join(feature_tables_dir, file)

            # Process the feature table for this subject
            process_subject(feature_table_qza, output_dir, subject_id)


if __name__ == "__main__":
    # Set up the argument parser for command-line usage
    parser = argparse.ArgumentParser(description="Calculate actual growth rates for multiple subjects.")
    parser.add_argument(
        "--feature_tables_dir", required=True, 
        help="Directory containing QIIME2 feature tables (.qza files) for all subjects."
    )
    parser.add_argument(
        "--output_dir", required=True, 
        help="Directory where outputs will be saved."
    )

    # Parse the arguments provided by the user
    args = parser.parse_args()

    # Run the main function with parsed arguments
    main(
        feature_tables_dir=args.feature_tables_dir,
        output_dir=args.output_dir
    )