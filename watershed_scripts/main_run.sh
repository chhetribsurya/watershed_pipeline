#!/bin/bash

# Define directories (if required)
# DATA_DIR="/scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep"

# Submit jobs for each piece of operation below:

# Calling expression outliers
sbatch piece1_job_array-final.sh

# Calling rare variants
sbatch piece2_job_array.sh

# Calling rare variant enrichmnent
sbatch piece3_job_array.sh
