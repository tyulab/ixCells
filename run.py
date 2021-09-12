#!/usr/bin/env python

from tasks import *
from utils.functions import get_folder
from pathlib import Path
import config

# Run all tasks here, can also comment out those not wanted
def main():
    # Create output.csv in output folder using files in data (starting with ASO)
    create_output() # get correct exp, zscore by bio replicate and mean of exp zscore (renormalization here)

    # # REQUIRED: output.csv
    create_avg_exp() # Create avg_exp.csv containing Exp zscore range and histograms from output.csv (unfiltered)
    create_control_hist() # Make control histograms

    # # REQUIRED: avg_exp_renormalized_to_neg.csv
    create_tiers(drop_output=False) # Create tiers.csv
    create_r_squared() # make r squared plots from WT, MT
    create_error_bars() # Make error bars

    # # REQUIRED: tiers.csv
    create_tier_plots() # Make plots showing tiers of samples


if __name__=="__main__":
    # Check all directories are created
    # renormalized
    folder = get_folder()
    Path(folder).mkdir(parents=True, exist_ok=True)
    Path(folder + "Control Plots").mkdir(parents=True, exist_ok=True)
    Path(folder + "Tiers").mkdir(parents=True, exist_ok=True)
    Path(folder + "Z Score Range Plots").mkdir(parents=True, exist_ok=True)
    # run tasks
    main()