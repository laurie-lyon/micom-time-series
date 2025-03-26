# Combine Simulated and Real Data for downstream analysis
# This script is to be used downstream of calculate_actual_growth_rates_*.py
# and simulated_growth_rates.py

#----------------------------------------------------------------------------
# Load in libraries
library(cowplot) # for plot grid
library(glue) # for string interpolation
library(lme4) # for linear modeling
library(tidyverse) # for data manipulation
library(ggh4x) # for expanded ggplot2 functionality

#----------------------------------------------------------------------------
# STEP 1: Read in all actual growth rate files and combine into 1 dataframe
#----------------------------------------------------------------------------
# OUTPUTS of calculate_actual_growth_rates_*.py are INPUT here in wide format
# create long df with all subjects included,
# add a column added for subject_id,
# and combine all subjects

#list of subject ids
subject_ids <- c("F01", "M01", "M02")

# Use map_dfr (from purrr package) to read, transform, and combine dataframes
actual_clr_long_combined <- map_dfr(subject_ids, function(subject) {
  # Read each subject's file (wide format dataframe) from folder
  df_wide <- read_csv(paste0("./data/actual_growth_rates_genus_clr/",
                             subject,
                             "_clr_actual_growth_rates_by_genus.csv"),
                      show_col_types = FALSE)
  # Transform to long format df and clean columns
  df_long <- df_wide %>%
    pivot_longer(cols = !Genus,
                 names_to = "sample_id",
                 values_to = "clr_change_abund") %>%
    rename(taxon = Genus) %>%
    #removes "g___" from taxon names for better readability
    mutate(taxon = str_sub(taxon, 4, -1),
           sample_id = as.numeric(sample_id),
           # add subject_id column directly
           subject_id = subject)
  return(df_long)
})

#----------------------------------------------------------------------------
# STEP 2: Read in all simulated growth rate files and combine into 1 dataframe
#----------------------------------------------------------------------------
# This section will read in all growth_rates.csv files from
#the unzipped folders and combine them into a single df (growth_df)

# OUTPUTS of simulated_growth_rates.py should be .zip files
#and unzipped folders with identical names
# these serve as INPUTS for this step

# grab list of all .zip files with the pattern we want
#(all should start with growth_)
zip_list <- Sys.glob("./data/growth_rates/growth*.zip")
# remove .zip from the end to get the unzipped folder basename
folder_list <- gsub(basename(zip_list), pattern = ".zip", replacement = "")

# establish pattern for file recognition in for loop
growth_file_pattern <- "./data/growth_rates/{folder_name}/growth_rates.csv"

# populate an empty list with dataframes
growth_df_list <- list()

for (folder_name in folder_list) {
  tmp_fp <- glue(growth_file_pattern)
  tmp_df <- read_csv(tmp_fp)
  tmp_df$folder_name <- folder_name
  # separate folder_name by "_" into new columns
  tmp_df <- separate(tmp_df,
    folder_name,
    # these column headers should correspond to folder_name labelling system
    into = c(
      "micom_step", # growth for "grow" output
      "subject_id",
      "model_db", # GSMN database (e.g. AGORA1, AGORA2, etc.)
      "solver", # e.g. osqp, gurobi, CPLEX
      "diet", # diet abbreviation specified in "diet_config.sh"
      "tradeoff" # coop tradeoff value set in simulate_growth_loop.sh
    ),
    sep = "_",
    remove = FALSE # if TRUE, removes original col
  )
  # create named list (each element named for its corresponding folder_name)
  # naming an element in a list uses [[double brackets]]
  growth_df_list[[folder_name]] <- tmp_df
}

# create a df that includes all data from all folders in folder_list
growth_df <- bind_rows(growth_df_list)

# remove "g___" from taxon names for better readability
micom_rates <- growth_df %>%
  mutate(taxon = str_sub(taxon, 4, -1)) %>%
  group_by(folder_name, sample_id)

#----------------------------------------------------------------------------
#(OPTIONAL) STEP 3: filter out only prevalent taxa
#----------------------------------------------------------------------------
# This step is optional and can be commented out if not desired
# This step is useful for reducing the number of taxa to plot

#set a prevalence threshold (0.50 = taxa present in 50% of samples)
prevalence_threshold <- 0.50

prevalent_micom_rates <- micom_rates %>%
  group_by(folder_name) %>%
  mutate(total_samples = n_distinct(sample_id)) %>%
  ungroup() %>%
  group_by(folder_name, taxon) %>%
  mutate(sample_counts = n_distinct(sample_id[!is.na(growth_rate)]),
         prevalence = sample_counts / total_samples) %>%
  ungroup() %>%
  filter(prevalence >= prevalence_threshold) %>%
  mutate(normed_epoch_time = (sample_id - min(sample_id)) /
           (max(sample_id) - min(sample_id))) %>%
  mutate(sample_id = as.numeric(sample_id))

#----------------------------------------------------------------------------
# STEP 4: Combine actual and simulated data
#----------------------------------------------------------------------------
# Combine actual and simulated data into a single dataframe (sim_real_data)
# Can use prevalent_micom_rates instead of micom_rates if desired

sim_real_data <- left_join(prevalent_micom_rates,
                           actual_clr_long_combined,
                           by = c("taxon", "sample_id", "subject_id")) %>%
  drop_na(growth_rate, clr_change_abund)

#----------------------------------------------------------------------------
# (OPTIONAL) STEP 5: Add metadata columns to sim_real_data
#----------------------------------------------------------------------------
# This step is optional and can be commented out if not desired
# Here, added sampling_date, sick_day, and date_vs_onset_illness columns
sim_real_data <- sim_real_data %>%
  #convert sample_id to POSIXct date format for cleaner plotting over time
  mutate(sampling_date = as.POSIXct(sample_id,
                                    origin = "1970-01-01",
                                    tz = "UTC")) %>%
  #add sick_day column for each subject (date of illness onset)
  mutate(sick_day = case_when(
                              subject_id == "F01" ~ 1202688000,
                              subject_id == "M01" ~ 1203465600,
                              subject_id == "M02" ~ 1203120000)) %>%
  # add column to differentiate samples taken before or after illness onset
  mutate(date_vs_onset_illness = case_when(
    sample_id < sick_day ~ "before",
    sample_id >= sick_day ~ "after"
  ))
#------------------------------------------------------------------------
# FINAL STEP: Save sim_real_data to a .csv file
#------------------------------------------------------------------------
# Save sim_real_data to a .csv file for downstream analysis
write_csv(sim_real_data, "./data/combined_sim_real_data_50.csv")
print("Combination complete. Data saved to ./data/combined_sim_real_data.csv")