import micom
import os
import pandas as pd
import zipfile
from pathlib import Path 
import argparse
from micom.workflows import build
from micom import Community
from micom.qiime_formats import load_qiime_medium
from micom.workflows import grow, save_results, complete_community_medium 

# Simulate growth rates for samples at each timepoint
# need to do this for each subject id


def load_subject_data(subject_id, qza_dir, collapse_on='genus'):
    """
    Load subject data and convert it to a MICOM-compatible format.
    Parameters:
    subject_id (str): The identifier for the subject.
    qza_dir (str): The directory where the QIIME2 artifact files are located.
    collapse_on (str, optional): The taxonomic level to collapse on. Default is 'genus'.
    Returns:
    micom.Community: A MICOM community object created from the subject's data.
    """
    
    feature_table_fp = os.path.join(qza_dir, f"{subject_id}_feature_table.qza")
    taxonomy_fp = os.path.join(qza_dir, f"{subject_id}_taxonomy.qza")
    subject_micom = micom.taxonomy.qiime_to_micom(feature_table_fp,
                                                  taxonomy_fp, 
                                                  collapse_on=collapse_on)
    
    return subject_micom

def add_suggested_metabolites(diet_og, diet_sugg, output_folder):
    """
    This function takes in the original diet and the micom suggested (completed) diet 
    and returns a new diet that includes the suggested metabolites 
    without removing the original ones.
    
    Inputs: 
    diet_og: pandas dataframe with the original diet
    diet_sugg: pandas dataframe with the diet from micom complete_community_medium
    output_folder: str, optional, directory where the added metabolites .csv file will be saved

    Returns:
    diet_new: pandas dataframe with the original and new nonzero elements of suggested diet
    """

    diet_og = diet_og.reset_index(drop=True)
    diet_sugg = diet_sugg.reset_index(drop=True)

    diet_merged = pd.merge(diet_og, diet_sugg, on=['reaction', 'metabolite'], how='outer', suffixes=('_og', '_sugg'))
    diet_merged["flux_diff"] = diet_merged["flux_sugg"] - diet_merged["flux_og"]
    added_metabolites = diet_merged[diet_merged["flux_diff"] > 0]
    added_metabolites = added_metabolites[["reaction", "metabolite", "global_id", "flux_sugg"]]
    added_metabolites = added_metabolites.rename(columns={"flux_sugg": "flux"})
    #write the added metabolites to a csv file 
    output_csv_fp = os.path.join(output_folder, "added_metabolites.csv")
    added_metabolites.to_csv(output_csv_fp, index=False)
    print(f"Added metabolites saved to {output_csv_fp}")
    #add added_metabolites to diet_og
    diet_new = pd.concat([diet_og, added_metabolites], ignore_index=True)
    #reindex diet_new
    diet_new = diet_new.reset_index(drop=True)
    return diet_new

def unzip_to_folder(growth_out_fp, out_folder):
    """
    Unzips the growth output file to the specified folder.
    Parameters:
    growth_out_fp (str): The path to the growth output file.
    out_folder (str): The folder to unzip the file to.
    """
    # unzip the growth output .zip file and save contents to a folder by the same name
    with zipfile.ZipFile(growth_out_fp, 'r') as zip_ref:
        zip_ref.extractall(out_folder)


def main(subject_id, qza_dir, 
         model_name, model_dir,
         pickled_gsmm_out, solver, 
         threads, diet_fp, 
         tradeoff, growth_out_fp):

    
    model_fp = os.path.join(model_dir, model_name)
    model_extract_fp = os.path.join(model_dir, Path(model_name).stem)

    subject_micom = load_subject_data(subject_id, qza_dir)

    diet_og = load_qiime_medium(diet_fp)
    #reindex diet_og to be row numbers [0:len(diet_og)]
    diet_og = diet_og.reset_index(drop=True)

    
    manifest = build(subject_micom,
                    out_folder=pickled_gsmm_out,
                    model_db=model_fp,
                    solver=solver,
                    threads=threads)
    
    diet_sugg = complete_community_medium(manifest, 
                                        model_folder=pickled_gsmm_out, 
                                        medium=diet_og, 
                                        community_growth=0.1, 
                                        min_growth=0.001, 
                                        minimize_components=True,
                                        max_import=1, 
                                        threads=threads)
    diet_sugg = diet_sugg.reset_index(drop=True)

    diet_new = add_suggested_metabolites(diet_og,
                                         diet_sugg,
                                         output_folder=pickled_gsmm_out)

    growth = grow(manifest, pickled_gsmm_out, 
                  medium=diet_new, tradeoff=tradeoff, 
                  threads=threads, presolve=True)
    save_results(growth, growth_out_fp)

    #unzip the growth output .zip file and save contents to a folder by the same name
    unzip_to_folder(growth_out_fp, growth_out_fp.replace(".zip", ""))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build and grow MICOM growth models")
    parser.add_argument("--subject_id", 
                        required=True, 
                        help="subject ID to process")
    parser.add_argument("--qza_dir",  
                        default="../data/qiime_outputs/", 
                        help="Path to .qza feature table")
    parser.add_argument("--model_dir", 
                        default="../data/models/", 
                        help="Path to model directory, also where .qza will be unzipped")
    parser.add_argument("--model_name", 
                        required=True, 
                        help="Name of .qza file for GSMM (e.g. agora103_genus.qza)")
    parser.add_argument("--pickled_gsmm_out", 
                        required=True, 
                        help="Output directory for GSMM .pickle files generated during build()")
    parser.add_argument("--solver", 
                        required=True,
                        default="osqp", 
                        help="Specify solver (e.g. osqp, gurobi, cplex)")
    parser.add_argument("--threads", 
                        type=int,
                        required=True, 
                        default=1,
                        help="Specify number of threads for paralellization ")
    parser.add_argument("--diet_fp", 
                        required=True, 
                        help="Path to qiime defined medium .qza (e.g. western diet gut agora)")
    parser.add_argument("--tradeoff", 
                        required=True, 
                        type=float,
                        help="Cooperative tradeoff (value between 0-1)")
    parser.add_argument("--growth_out_fp", 
                        required=True, 
                        help="Path for output growth.zip from micom grow()")
    

    args = parser.parse_args()

    main(args.subject_id, args.qza_dir, 
        args.model_name, args.model_dir,
        args.pickled_gsmm_out, args.solver, 
        args.threads, args.diet_fp, 
        args.tradeoff, args.growth_out_fp)

