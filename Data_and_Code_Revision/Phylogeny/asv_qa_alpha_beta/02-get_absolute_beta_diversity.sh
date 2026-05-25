#!/bin/bash
set -e pipefail

# ******************************************************************************
# File: a.sh
# Author: Mingxing Wang
# Email: xing592798030@163.com
# Date: 2025-12-25 15:09:52
# License: Copyright (C) 2025 Mingxing Wang. All rights reserved.
# Reference: Mingxing Wang
# Description: 
# ******************************************************************************


# Import configuration
source "$config/setup.sh"


# Define variable
env=r451


# Define function


# Execute main process
conda activate "$env"


export OPENBLAS_NUM_THREADS=1

for kin in 16S ITS Protist
do
   	Rscript --vanilla 02-get_absolute_beta_diversity.R $kin >> log/02-get_absolute_beta_diversity.log 2>&1 &
done    



