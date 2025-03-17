import os
import re

# Define your directory path
dir_path = "../data/growth_rates/"

# Define regex pattern:
# Matches files starting with "growth_" that don't have agora201 as the c parameter
pattern = re.compile(r"^(growth_[^_]+_)(?!agora201_)([^_]+_[^_]+_[^_]+)$")

# Iterate over the files
for filename in os.listdir(dir_path):
    # Check for regex match
    match = pattern.match(filename)
    if match:
        prefix = match.group(1)
        suffix = match.group(2)
        
        # Construct new filename inserting agora1 as the missing parameter c
        new_filename = f"{prefix}agora1_{suffix}"
        
        # Rename the file
        old_path = os.path.join(dir_path, filename)
        new_path = os.path.join(dir_path, new_filename)
        
        print(f"Renaming: {filename} -> {new_filename}")
        os.rename(old_path, new_path)
