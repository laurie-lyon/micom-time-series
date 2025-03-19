# Combine Simulated and Real Data for downstream analysis
# This script is to be used downstream of calculate_actual_growth_rates_*.py
# and simulated_growth_rates.py

# Load in libraries
library(cowplot) # for plot grid
library(glue) # for string interpolation
library(lme4) # for linear modeling
library(tidyverse) # for data manipulation
library(ggh4x) # for expanded ggplot2 functionality

# OUTPUTS of simulated_growth_rates.py should be .zip files
#and unzipped folders with identical names
# This section will read in all growth_rates.csv files from
#the unzipped folders and combine them into a single dataframe

# grab list of all .zip files with the pattern we want
#(all should start with growth_)
zip_list <- Sys.glob("../data/growth_rates/growth*.zip")
# remove .zip from the end to get the unzipped folder basename
folder_list <- gsub(basename(zip_list), pattern = ".zip", replacement = "")

# establish pattern for file recognition in for loop
growth_file_pattern <- "../data/growth_rates/{folder_name}/growth_rates.csv"

# populate an empty list with dataframes
growth_df_list <- list()

for (folder_name in folder_list) {
    tmp_fp <- glue(growth_file_pattern)
    tmp_df <- read_csv(tmp_fp)
    tmp_df$folder_name <- folder_name
    # separate folder_name by "_" into new columns
    tmp_df <- separate(tmp_df,
        folder_name,
        into = c(
            "micom_step",
            "subject_id",
            "model_db",
            "solver",
            "diet",
            "tradeoff"
        ),
        sep = "_",
        remove = FALSE # if TRUE, removes original col
    )
    # create named list (each element named for its corresponding folder_name)
    # naming an element in a list uses [[double brackets]]
    growth_df_list[[folder_name]] <- tmp_df
    # print(head(tmp_df))
}

# create a df that includes all data from all folders in folder_list
growth_df <- bind_rows(growth_df_list)


