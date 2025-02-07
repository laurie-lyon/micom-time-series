"""
time_series_data_wrangling_w_parser.py

This script processes OTU tables and metadata for microbiome analysis by:
1. Filtering OTU data and metadata for specific subjects (anonymized names).
2. Averaging microbial abundances for each OTU across samples if mulpiple samples are present for a single day.
3. Identifying missing daily timepoints and interpolating OTU abundance by Piecewise Cubic Hermite Interpolating Polynomial (PCHIP).
4. Combining interpolated time points with real/averaged daily timepoints and preparing them for QIIME2.
5. Updating metadata to include interpolated samples by generating sample IDs and data types ("Real" or "Interpolated").

The script uses command-line arguments for file paths, output directories, and subject IDs, making it flexible for different datasets.

**Outputs**:
- Combined OTU tables for QIIME2 in `.tsv` format.
- Updated metadata files with real and interpolated sample entries labeled accordingly.

Authors: Laurie Lyon and Casey Martin
Date: 01/21/2025
"""

# Import required libraries
from biom import load_table
import pandas as pd
import numpy as np
from scipy.interpolate import PchipInterpolator
from pathlib import Path
import argparse

# ---- Helper Functions ----

def biom2df(biom_table):
    """
    Converts a BIOM table into a pandas DataFrame.

    Parameters:
    - biom_table: A BIOM-format table containing OTU labels, relative abundances, and sample ids 
    
    Output:
    - a pandas datafram with OTU labels as rows and sample ids as columns
    """
    data_matrix = biom_table.matrix_data.toarray()  # Convert sparse matrix to dense format
    df = pd.DataFrame(data_matrix,
                      index=biom_table.ids(axis="observation"),
                      columns=biom_table.ids(axis="sample")) # Create pandas dataframe 
    return df


def meta_biom_filter(metadata, biom_df, anon_name, sample_col="#SampleID", subject_col="ANONYMIZED_NAME"):  
    """
    Filters metadata and OTU table for a specific subject.

    Parameters:
    - metadata: Metadata as a pandas DataFrame.
    - biom_df: OTU abundance DataFrame converted from a BIOM table.
    - anon_name: Subject ID to filter by.
    - sample_col: Column name for sample IDs in metadata, default is "#SampleID".
    - subject_col: Column name for subject IDs in metadata, default is "ANONYMIZED_NAME".

    Returns:
    - Filtered metadata (subject_metadata) and OTU table (subject_biom) for the specified subject.

    These outputs are used as inputs for the downstream avg_sample_by_day and interp_missing_day functions.

    """
    # Set index to sample ID
    metadata.index = metadata[sample_col]
    # subset metadata for subject id
    subject_metadata = metadata[metadata[subject_col] == anon_name]

    # Biom table is columns: samples, rows: OTUs
    # We need to get the sample ids (columns) from the biom table that are also 
    # in the metadata under the sample_col (usually #SampleID) header.
    #
    # We also need to get only the samples that are in the metadata that are also present
    # in the biom table.
    # This is the intersection of the set of samples in the metadata and the set of samples in the biom table.
    # https://en.wikipedia.org/wiki/Set_(mathematics)#Basic_operations
    shared_indices = set(subject_metadata[sample_col]) & set(biom_df.columns) # set intersection
    shared_indices = list(shared_indices) # convert back to list because we cannot index dataframes using sets

    # filter subject-specific biom tables and metadata for only shared samples (shared_indices)
    subject_biom = biom_df[shared_indices]
    subject_metadata = subject_metadata.loc[shared_indices]

    # return filtered metadata and biom table specific for each subject with only shared indices present
    return subject_metadata, subject_biom


def avg_sample_by_day(subject_metadata, 
                      subject_biom, 
                      sample_col="#SampleID", 
                      time_col="epoch_time", 
                      subject_col="ANONYMIZED_NAME"):
    """
    Averages OTU abundances per day for each subject.

    Parameters: 
    - subject_metadata: Metadata for a specific subject (created by meta_biom_filter).
    - subject_biom: OTU table for a specific subject (created by meta_biom_filter).
    - sample_col: Column name for sample IDs in metadata, default is "#SampleID".
    - time_col: Column name for time points in metadata, default is "epoch_time".
    - subject_col: Column name for subject IDs in metadata, default is "ANONYMIZED_NAME".

    Returns:
    - avg_subject_metadata: Metadata for each subject with one entry per day (if that day has samples).
    - avg_day: OTU table with average abundances per day for each subject (if that day has samples)

    Outputs from this function (avg_subject_metadata, avg_day) are used as inputs for the interp_missing_day function.

    """
    # convert the time column in subject_metadata to integers for sorting
    subject_metadata[time_col] = subject_metadata[time_col].astype(int)

    # get unique timepoints for each subject and sort by time
    timepoints = sorted(subject_metadata[time_col].unique())

    # calculate average abundance for each day (avg_day)
    avg_day = [subject_biom[subject_metadata.loc[subject_metadata[time_col] == day, sample_col]].mean(axis=1) for day in timepoints]

    #index the avg_day dataframe by the unique timepoints for each subject
    avg_day = pd.DataFrame(avg_day, index=timepoints)

    # indexes the subject_metadata by the unique timepoints for each subject
    avg_subject_metadata = subject_metadata.groupby(time_col).first().reset_index()
    
    return avg_subject_metadata, avg_day


def epoch_datetime_range(timepoints): 
    """
    Generates a range of epoch times (one per day) between the earliest and latest timepoints.

    Output: 
    - epoch_time_range: A range of epoch times (in seconds) between the earliest and latest timepoints

    This function is dependent on the avg_sample_by_day function where timepoints is defined.
    epoch_time_range is used as an input for the missing_timepoints function.

    """
    # finds the first and last timepoints in the timepoints list
    earliest = pd.to_datetime(min(timepoints), unit='s')
    latest = pd.to_datetime(max(timepoints), unit='s')

    # creates a range of dates (by day) between the earliest and latest timepoints
    date_range = pd.date_range(start=earliest, end=latest, freq='D')

    # converts the date range to integers in pandas and converts to seconds for epoch time
    epoch_time_range = date_range.astype('int64') // 10**9  # Convert to seconds

    return epoch_time_range


def missing_timepoints(epoch_time_range, timepoints):
    """
    Identifies timepoints missing from the dataset.

    - output: a list of missing timepoints

    This function is dependent on the epoch_datetime_range function.
    The output list is used as an input for the interp_missing_day function.

    """
    return list(set(epoch_time_range) - set(timepoints))


def interp_missing_day(avg_subject_meta, avg_subject_biom):
    """
    Interpolates missing daily time points in an OTU table.

    Inputs: 
    - avg_subject_meta: Metadata for each subject with one entry per day (generated by avg_sample_by_day function).
    - avg_subject_biom: OTU table with average abundances per day for each subject (generated by avg_sample_by_day function).
    
    """

    pchip = PchipInterpolator(avg_subject_biom.index, avg_subject_biom)
    subject_date_range = epoch_datetime_range(avg_subject_meta["epoch_time"])
    subject_missing_times = missing_timepoints(subject_date_range, avg_subject_meta["epoch_time"])
    interp_biom = pchip(subject_missing_times)
    interp_biom = pd.DataFrame(interp_biom, index=subject_missing_times, columns=avg_subject_biom.columns)

    return interp_biom


# ---- Main Script ----
def main(input_biom_path, metadata_path, output_dir, final_output_dir, subject_ids):
    """
    Main pipeline for processing OTU tables and metadata for microbiome analysis.
    """
    # check if directory exists and make a directory if one doesn't already exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(final_output_dir).mkdir(parents=True, exist_ok=True)

    print("Loading BIOM table and metadata...")
    myotubiom = load_table(input_biom_path)
    mymeta = pd.read_csv(metadata_path)
    mybiom = biom2df(myotubiom)

    output_interp_bioms = {}
    output_avg_subject_meta = {}
    output_avg_subject_biom = {}

    for subject in subject_ids:
        print(f"Processing data for subject {subject}...")
        subject_metadata, subject_biom = meta_biom_filter(mymeta, mybiom, subject)
        avg_subject_meta, avg_subject_biom = avg_sample_by_day(subject_metadata, subject_biom)
        subject_interp_biom = interp_missing_day(avg_subject_meta, avg_subject_biom)

        subject_interp_biom.to_csv(f"{output_dir}{subject}_interp_biom.csv")
        avg_subject_biom.to_csv(f"{output_dir}{subject}_avg_biom.csv")
        avg_subject_meta.to_csv(f"{output_dir}{subject}_avg_meta.csv")

        output_interp_bioms[subject] = subject_interp_biom
        output_avg_subject_meta[subject] = avg_subject_meta
        output_avg_subject_biom[subject] = avg_subject_biom

    for subject in subject_ids:
        print(f"Combining and transposing OTU tables for subject {subject}...")
        combined_otu = pd.concat([output_avg_subject_biom[subject], output_interp_bioms[subject]], axis=0).transpose()
        combined_otu.index.name = "#OTU ID"
        combined_otu.to_csv(f"{final_output_dir}{subject}_combined_otu_qiime2.tsv", sep="\t")

    for subject in subject_ids:
        print(f"Updating metadata for subject {subject}...")
        interp_metadata = pd.DataFrame({
            "#SampleID": "1015_" + pd.to_datetime(output_interp_bioms[subject].index, unit="s").strftime("%m_%d") + f"_{subject}_interp",
            "epoch_time": output_interp_bioms[subject].index,
            "ANONYMIZED_NAME": subject,
            "data_type": "Interpolated"
        })
        output_avg_subject_meta[subject]["data_type"] = "Real"
        updated_metadata = pd.concat([output_avg_subject_meta[subject], interp_metadata]).sort_values(by="epoch_time")
        updated_metadata.to_csv(f"{final_output_dir}{subject}_updated_metadata.csv", index=False)


# ---- Command-Line Interface ----
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process OTU tables and metadata for microbiome analysis.")
    parser.add_argument("-i", "--input_biom", required=True, help="Path to input BIOM file.")
    parser.add_argument("-m", "--metadata", required=True, help="Path to metadata file.")
    parser.add_argument("-n", "--intermediate_output_dir", required=True, help="Directory for intermediate outputs.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory for final outputs.")
    parser.add_argument("-l", "--subject_ids", required=True, nargs="+", help="List of subject IDs to process.")
    args = parser.parse_args()

    main(args.input_biom, args.metadata, args.intermediate_output_dir, args.output_dir, args.subject_ids)