#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
from tqdm import tqdm

from utils.functions import *
from utils.plots import *
from utils.stats import *

# limit on z score range to filter out before making list of tiers
TIERS_THRESHOLD = 1.5
# threshold for output
Z_SCORE_THRESHOLD = 1.8

# make output file
# steps- renormalize ddCt, recalculate Exp, calculate abs z scores, drop outliers, avg z scores on samples, assign plates, avg on plates
def create_output():
    # specify path to folder containing all csvs and plate sheet + output to be ignored
    path = os.getcwd()+"\data"
    plate_file = "ixCells_Round 1_2021-06-22_TN09_551ASOs_plate id adjusted.csv"
    output_file = "output.csv"
    csv = get_csv(path, plate_file)

    pd.set_option('display.max_columns', None)

    df_list = []
    print('found ' + str(len(csv)) + ' csv files')
    # read all csvs into dataframes and concatenate
    for f in tqdm(csv):
        if plate_file not in f:
            df = pd.read_csv(f, encoding='latin-1')
            df.columns = df.columns.str.strip()
            # rename short description to microsynth id to have consistent columns, original naive normalized cols
            df = df.rename(columns={'ASO Short Description': 'ASO Microsynth ID'})
            # drop fully empty rows
            df = df.dropna(how='all').reset_index(drop=True)
            # redo dCt from Crossing Point
            df = redo_dct(df)
            # drop irrelevant columns
            df = df[['Experiment Name', 'Position', 'SampleName', 'ASO Microsynth ID', 'dCt', 'ddCt', 'Exp']]

            # renormalizing
            df = df.rename(columns={'ddCt': 'ddCt from Naive', 'Exp':'Exp from Naive'})
            neg_ddct(df)


            df = df.replace(to_replace='^Na.*e', value='Naive', regex=True)
            df_list.append(df)

    df = pd.concat(df_list, ignore_index=True)
    # pd.set_option('display.max_rows', df.shape[0]+1)

    # Calculate exp
    calc_exp(df)
    # renormalizing
    calc_exp(df, col='ddCt from Naive', output_col='Exp from Naive')

    # Calculate z scores for Exp
    abs_zscore(df)
    # drop outliers
    # set threshold on z exp zscores
    threshold = flag_outliers(df, col='Exp_zscore', greater_than=Z_SCORE_THRESHOLD)
    dropped = df[threshold].reset_index(drop=True)
    dropped.to_csv("output/output_dropped.csv")
    df[~threshold].reset_index(drop=True).to_csv("output/output_filtered.csv")
    df.loc[threshold,'Exp'] = np.nan
    # df = df[~threshold].reset_index(drop=True)
    # Average z scores
    # todo: remove
    # avg_zscore(df)
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

    # calculate avg exp, range
    df = exp_zscore_range(df)
    df.to_csv("output/avg_exp_renormalized_to_neg.csv", index=False)
    # make histogram from ranges
    type_hist(df, 'Avg Exp_zscore range')

def create_tiers():
    # default path for table to read from
    table_file = "output/avg_exp_renormalized_to_neg.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # remove 3mmol and total from view
    df = remove_3(df)
    df = remove_total(df)

    # flag ranges > threshold and drop
    df = df.dropna(subset=['Avg Exp_zscore range']).reset_index(drop=True)
    ranges = flag_outliers(df, col='Avg Exp_zscore range', greater_than=TIERS_THRESHOLD)
    dropped = df[ranges].sort_values('Avg Exp_zscore range', ascending=True).reset_index(drop=True)
    dropped_groups = dropped.groupby('SampleName').first().sort_values(['SampleName'], na_position='last')
    dropped_groups = dropped_groups[['ASO Microsynth ID', 'Test plate #']]
    dropped_groups['Threshold: ' + str(TIERS_THRESHOLD)] = np.nan
    dropped_groups.to_csv("output/tiers_dropped.csv", index=True)

    # remove anything in dropped samples from tiering
    dropped_samples = dropped['SampleName'].unique()
    df = df[~df['SampleName'].isin(dropped_samples)].reset_index(drop=True)

    tiers = tierlist(df)

    tiers.sort_values(['Tier','SampleName'], ascending=True, inplace=True)
    tiers.to_csv("output/tiers.csv", index=True)


# make control histograms
def create_control_hist():
    # default path for table to read from
    table_file = "output/output.csv"
    plate_output = "plates_table.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # TODO: flag any controls across plates > 1.55
    print('control plots:')
    outliers = flag_outliers(df, col='Exp', greater_than=1.55)
    df = df[~outliers].reset_index(drop=True)
    # histogram
    control_hist(df)

    # # drop columns/nan rows and export
    # plates = drop_nan(df)
    # plates.to_csv("output/" + plate_output, index=False)


def create_tier_plots():
    file = "output/tiers.csv"
    # lst of column names which needs to be string
    lst_str_cols = ['Tier']
    # use dictionary comprehension to make dict of dtypes
    dict_dtypes = {x: 'str' for x in lst_str_cols}
    # use dict on dtypes
    df = pd.read_csv(file, dtype=dict_dtypes)
    # scatter plot
    tier_scatter(df)
    # histogram
    # tier_hist(df)

# r squared analysis
def run_r_squared():
    file = "output/avg_exp_renormalized_to_neg.csv"
    # store in df
    df = pd.read_csv(file, encoding='latin-1')
    pd.set_option('display.max_columns', None)
    # keep only 10mmol
    df = df[df['SampleName'].str.contains('_10')]
    # get replicates
    for type in ['WT','MT']:
        type_df = get_type(df, type=type)
        b1,b2 = get_replicates(type_df)
        r2_plot(b1,b2,title=type)

def create_box_plot():
    file = "output/avg_exp_renormalized_to_neg.csv"
    # store in df
    df = pd.read_csv(file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # test versions
    # box_plot(df)
    error_plot(df)
    # lm_plot(df)
    # bar_plot(df)