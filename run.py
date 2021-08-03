#!/usr/bin/env python

from tasks import *

# Run all tasks here, comment out those not wanted
def main():
    # Create output.csv in output folder using files in data (starting with ASO)
    # Get renormalized output from data, exp zscore and avg exp zscore
    # create_output()

    # # REQUIRED: output.csv
    # Create avg_exp_renormalized_to_neg.csv containing Exp zscore range and histograms from output.csv (unfiltered)
    # create_avg_exp()
    # Make control histograms
    # create_control_hist()

    # # REQUIRED: avg_exp_renormalized_to_neg.csv
    # Create tiers.csv
    # create_tiers()
    # make r squared plots from WT, MT
    # run_r_squared()
    # Make box plots
    create_box_plot()

    # # REQUIRED: tiers.csv
    # Make plots showing tiers of samples
    # create_tier_plots()



if __name__=="__main__":
    main()