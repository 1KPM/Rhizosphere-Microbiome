#!/bin/bash
set -e pipefail

# ******************************************************************************
# File: 03-get_top10order_network_property_linux.sh
# Author: Mingxing Wang
# Email: xing592798030@163.com
# Date: 2026-03-05 17:05:59
# License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
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

for typ in all inter intra
do
	Rscript --vanilla 03-robustness_linux.R $typ >> log/03-robustness_linux_${typ}.log 2>&1 &
done

