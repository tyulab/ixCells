#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import dose_response
import glob
import os

from functions import get_csv, flag_outliers, assign_plates
from histograms import type_hist, control_hist
from stats import neg_ddct, calc_exp, abs_zscore, avg_zscore, avg_plates, exp_zscore_range, tierlist

# make output file
# steps- renormalize ddCt, recalculate Exp, calculate abs z scores, drop outliers, avg z scores on samples, assign plates, avg on plates
def create_output():
    # TODO: parser?
    # specify path to folder containing all csvs and plate sheet + output to be ignored
    path = os.getcwd()+"\data"
    plate_file = "ixCells_Round 1_2021-06-22_TN09_551ASOs_plate id adjusted.csv"
    output_file = "output.csv"
    csv = get_csv(path, plate_file)

    pd.set_option('display.max_columns', None)

    df_list = []
    print('found ' + str(len(csv)) + ' csv files')
    # read all csvs into dataframes and concatenate
    for f in csv:
        if plate_file not in f:
            df = pd.read_csv(f, encoding='latin-1')
            df.columns = df.columns.str.strip()
            # rename short description to microsynth id to have consistent columns, original naive normalized cols
            df = df.rename(columns={'ASO Short Description': 'ASO Microsynth ID'})
            # drop irrelevant columns
            df = df[['Experiment Name', 'Position', 'SampleName', 'ASO Microsynth ID', 'dCt', 'ddCt', 'Exp']]
            df = df.rename(columns={'ddCt': 'ddCt from Naive', 'Exp':'Exp from Naive'})
            # calculate ddCt
            neg_ddct(df)
            df = df.replace(to_replace='^Na.*ve', value='Naive', regex=True)
            df_list.append(df)

    df = pd.concat(df_list, ignore_index=True)
    # pd.set_option('display.max_rows', df.shape[0]+1)
    # Calculate exp
    calc_exp(df)
    calc_exp(df, col='ddCt from Naive', output_col='Exp from Naive')
    # Calculate z scores for Exp
    abs_zscore(df)
    # drop outliers
    # set threshold on z exp zscores
    threshold = flag_outliers(df, col='Exp_zscore', greater_than=1.8)
    dropped = df[threshold].reset_index(drop=True)
    dropped.to_csv("output/output_dropped.csv")
    # Average z scores
    avg_zscore(df)
    assign_plates(df, plate_file)
    avg_plates(df) # can separate this part later?
    print(df)
    # export to output file
    df.to_csv("output/"+output_file, index=False)

# create avg exp file
def create_avg_exp():
    # default path for table to read from
    table_file = "output/output.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    df = df[['Experiment Name','SampleName','ASO Microsynth ID','Exp_zscore','Exp','Test plate #']]
    pd.set_option('display.max_columns', None)

    # functions
    df = exp_zscore_range(df)
    df.to_csv("output/avg_exp_renormalized_to_neg.csv", index=False)

    type_hist(df, 'Avg Exp_zscore range')

def create_tiers():
    # default path for table to read from
    table_file = "output/avg_exp_renormalized_to_neg.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # flag >1.55 zscore range and/or drop
    df = df.dropna(subset=['Avg Exp_zscore range']).reset_index(drop=True)
    threshold = 1.25
    ranges = flag_outliers(df, col='Avg Exp_zscore range', greater_than=threshold)
    dropped = df[ranges].sort_values('Avg Exp_zscore range', ascending=True).reset_index(drop=True)
    dropped_samples = dropped.groupby('SampleName').first().sort_values(['SampleName'], na_position='last')
    dropped_samples = dropped_samples[['Experiment Name', 'ASO Microsynth ID', 'Test plate #', 'Avg Exp_zscore range']]
    dropped_samples['Threshold: ' + str(threshold)] = np.nan
    dropped_samples.to_csv("output/tiers_dropped.csv", index=False)
    # remove anything weihtout range information
    df = df[~ranges].reset_index(drop=True)
    df = tierlist(df)
    grouped_samples = df.groupby('SampleName').first().sort_values(['Tier', 'SampleName'], na_position='last')
    grouped_samples = grouped_samples[['ASO Microsynth ID', 'Test plate #', 'Tier']]
    grouped_samples.to_csv("output/tiers.csv", index=False)
    # drop ranges


# make control histograms
def create_control_hist():
    # default path for table to read from
    table_file = "output/output.csv"
    plate_output = "plates_table.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # TODO: flag any controls across plates > 1.55
    outliers = dose_response.flag_outliers(df, col='Exp', greater_than=1.55)
    df = df[~outliers].reset_index(drop=True)
    # histogram
    control_hist(df)

    # # drop columns/nan rows and export
    # plates = drop_nan(df)
    # plates.to_csv("output/" + plate_output, index=False)