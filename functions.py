#!/usr/bin/env python
import re

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# count outliers and return
def flag_outliers(df, col='Exp_zscore', greater_than=2.0):
    count = df[df[col] > greater_than].count()[0]
    # df[zscore_col+'<'+greater_than] = True
    # df.loc[df[col] > greater_than, exp_col] = False
    print("%s: %d greater than %s" % (col, count, greater_than))
    return df[col] > greater_than


def clear_outliers(df, zscore_col='Exp_zscore', exp_col='Exp', greater_than=1.8):
    count = df[df[zscore_col] > greater_than].count()[0]
    df.loc[df[zscore_col] > greater_than, exp_col] = np.nan
    print("%s: %d greater than %s" % (zscore_col, count, greater_than))
    return df

# get list of all csvs
def get_csv(path, plate_file):
    # get all csvs in path holding data, ignore files that don't start with ASO
    csv = [f for f in glob.iglob(os.path.join(path, "**/*.csv"),recursive=True) \
           if os.path.basename(f).startswith('ASO')]
    return csv

# read plates and add column to output
def assign_plates(df, file="ixCells_Round 1_2021-06-22_TN09_551ASOs_plate id adjusted.csv"):
    plates = pd.read_csv(file, encoding='latin-1')
    plates.columns = plates.columns.str.strip()
    plates = plates[['Test plate #', 'ASO Microsynth ID']]
    df["Test plate #"] = np.nan
    for i in range(len(df)): # try to find _P[0-9]+
        # sample = str(df.loc[i]['SampleName']).split("_")
        find_plate_no = re.search('_P\d+', str(df.loc[i,'SampleName'])) # search for plate number on sample name
        if find_plate_no is None:
            find_plate_no = re.search('_P\d+', str(df.loc[i,'ASO Microsynth ID']))  # search for plate number on microsynth id
        # add plate #s to naive and control
        if find_plate_no is not None:
        # if "Ionis" in sample[0] or "Naive" in sample[0]:  # check if control or naive
            # get plate number from end eg _P##
            plate_no = find_plate_no.group(0)[2:]
            df.loc[i, 'Test plate #'] = int(plate_no)
        else:
            # else attempt to find using table
            aso = str(df.loc[i]['ASO Microsynth ID']).split("_")[0]
            if aso in plates['ASO Microsynth ID'].values:
                df.loc[i,'Test plate #'] = int(plates[plates['ASO Microsynth ID'] == aso]['Test plate #'])


def hide_naive_control(df):
    df = df[~df['SampleName'].str.contains('Naive|Ionis', case=False)].copy().reset_index(drop=True)
    return df


def get_avg_exp(df, sample_name):
    try:
        avg_exp = df[sample_name]
    except:
        avg_exp = None
    return avg_exp