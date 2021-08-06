#!/usr/bin/env python
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy import stats

# Functions for commonly used math and all calculations

# get absolute value z score from dataframe on column
from utils.functions import hide_naive_control, get_avg_exp
import config


def abs_zscore(df, col='Exp'):
    # z score lambda function
    zscore = lambda x: abs((x - x.mean()) / x.std(ddof=1))
    # apply to exp column based on sample name]
    # TODO: remove std (for testing)
    df[col +'_std'] = df[col].groupby(df.index // 3).transform('std')
    df[col+'_zscore'] = df.groupby(['Experiment Name','SampleName'])[col].transform(zscore)

# redo dCt based on Crossing Point
def redo_dct(df, CP_1='CrossingPoint', CP_2='CrossingPoint.1'):
    df['dCt'] = abs(df[CP_1]-df[CP_2])
    return df

# redo ddCt normalized on neg control
# TODO: make it to apply after merging sheets, low priority
def neg_ddct(df):
    # get all values Ionis676.., group by experiment name/sample name
    ionis_neg = df[df['SampleName'].str.contains('(?:Ionis676630|Ionis 676630).*_10', na=False)]
    # just more convenient to keep all columns
    df['Avg dCt Ionis676630'] = ionis_neg['dCt'].mean()
    # do difference
    df['ddCt'] = df['dCt'] - df['Avg dCt Ionis676630']

# redo exp based on ddct column
def calc_exp(df, col='ddCt', output_col='Exp'):
    expression = lambda x: 2 ** (-x)
    df[output_col] = df[col].transform(expression)

# use unbiased estimator
def calc_std(df, col='Exp', output_col='Exp std'):
    df['Exp_std'] = df.groupby(['Experiment Name', 'SampleName'])['Exp'].transform('std')
    return df

def calc_sem(df, col='Exp', output_col='Exp SEM'):
    SEM = lambda x: x.std() / math.sqrt(x.count())
    df['Exp_SEM'] = df.groupby(['Experiment Name', 'SampleName'])['Exp'].transform(SEM)
    return df


# avg zscore by plate
def avg_plates(df):
    # fixed: includes naive and control into calculations
    # sort on test plate and sample
    # df["Test plate #"] = pd.to_numeric(df["Test plate #"], errors='coerce')
    df.sort_values(['Test plate #', 'SampleName'],ascending=True, na_position='last', ignore_index=True, inplace=True)
    # get average zscore by plate value
    plates = df['Exp_zscore'].groupby(df['Test plate #'])
    df["Avg_plate"] = plates.transform('mean')
    # get min max range based on exp_zscore values
    df["Min_plate"] = plates.transform("min")
    df["Max_plate"] = plates.transform("max")
    # set duplicate values to nan
    duplicate = df.duplicated(['Test plate #', 'Avg_plate'],keep='first')
    df.loc[duplicate, ['Avg_plate', "Min_plate", "Max_plate"]] = np.nan

# avg zscore by bio replicates
def avg_zscore(df, col='Exp_zscore'):
    # average by replicate
    df["Avg "+col] = np.nan
    # print(df[col].groupby(df.index // 3).mean())
    mean = df[col].groupby(df.index // 3).transform('mean')
    df.loc[::3,'Avg '+col] = mean

    # change naive
    df["Avg_naive"] = np.nan
    # transform mean on naive groups
    df.loc[df['SampleName'].str.contains("Naive", case=False),"Avg_naive"] =\
        df[df['SampleName'].str.contains("Naive", case=False)].groupby(['Experiment Name','SampleName'])['Exp_zscore'].transform('mean')
    # set duplicate values to nan -- ignore missing microsynth id
    duplicate = df.duplicated(['Experiment Name','SampleName', 'Avg_naive'], keep='first')
    df.loc[df['SampleName'].str.contains("Naive", case=False) & duplicate, 'Avg_naive'] = np.nan

# pass in column, mean, and std of same row length (use transform)
def zscore(df, col, mean, std):
    df[col+'_zscore'] = df[col].sub(mean).div(std)# .abs()

# ttest
def two_sample_ttest(df, col):
    df[col+' t_stat'] = df[col + ' p_value'] = df[col + ' dof'] = np.nan
    # iterate over every 6 rows
    for i in range(0, len(df), 6):
        t, p, dof = welch_ttest(df[i:i+3][col], df[i+3:i+6][col])
        df.loc[i,[col+' t_stat', col + ' p_value', col + ' dof']] = [t, p, dof]
    # print(df)

# helper fn for t test
def welch_ttest(x, y):
    ## Welch-Satterthwaite Degrees of Freedom ##
    t, p = stats.ttest_ind(x, y, equal_var=False, nan_policy='omit')
    return t, p, welch_dof(x, y)

# helper fn for t test
def welch_dof(x, y):
    dof = (x.var() / x.size + y.var() / y.size) ** 2 / ((x.var() / x.size) ** 2 / (x.size - 1) + (y.var() / y.size) ** 2 / (y.size - 1))
    return dof

# calculate average zscore by replicate and range with output (unfiltered)
def exp_zscore_range(df):
    # average by replicate
    df["Avg Exp"] = np.nan
    # print(df['Exp_zscore'].groupby(df.index // 3).mean())
    mean = df['Exp'].groupby(df.index // 3).transform('mean')
    df.loc[::3, 'Avg Exp'] = mean
    # get mean and std over each sample (every 6 replicates) to use on Avg Exp
    mean = df.groupby(['Experiment Name','SampleName'])['Exp'].transform('mean')
    std = df.groupby(['Experiment Name','SampleName'])['Exp'].transform('std')
    # zscore on Avg Exp
    zscore(df, 'Avg Exp', mean, std)

    # get range of each group
    range = lambda x: abs(x.max() - x.min())
    # apply to exp column based on sample name
    # print(df.groupby(['Experiment Name','SampleName']).groups.keys())
    df['Avg Exp_zscore range'] = df.groupby(['Experiment Name', 'SampleName'])['Avg Exp_zscore'].transform(range)

    # two_sample_ttest(df, 'Exp') # confirm results
    return df

# helper fn: pass MT and WT, assign to tiers
def create_tier(df, col='Avg Exp', new_col='Tier'):
    df[new_col] = np.vectorize(config.DICT.get)(np.digitize(df[col], config.BINS, right=False))


def tierlist(df):
    # get mean expression of sample
    df['Avg Exp'] = df.groupby(['Experiment Name','SampleName'])['Exp'].transform('mean')
    # create MT and WT views
    df_mt = df[df['Experiment Name'].str.contains('MT', case=False)].reset_index(drop=True)[['SampleName','ASO Microsynth ID','Test plate #','Avg Exp']]
    df_wt = df[df['Experiment Name'].str.contains('WT', case=False)].reset_index(drop=True)[['SampleName','ASO Microsynth ID','Test plate #','Avg Exp']]
    df = df[['SampleName','ASO Microsynth ID','Avg Exp']]
    # df['Tier'] = np.nan
    # group by name and keep only first row #todo: calculate mean here instead?
    df_mt = df_mt.groupby(df_mt['SampleName']).first()
    # rename avg exp in MT
    df_mt = df_mt.rename(columns={'Avg Exp': 'MT Avg Exp'})
    df_wt = df_wt.groupby(df_wt['SampleName']).first()
    # save views
    # df_mt.to_csv(get_folder()+"/df_mt.csv")
    # df_wt.to_csv(get_folder()+"/df_wt.csv")
    # TODO: make sure count of both is correct
    # add WT Avg exp column to df_mt
    df_mt['WT Avg Exp'] = df_wt['Avg Exp']

    # create tiers for both df_mt, df_wt views
    create_tier(df_mt, col='MT Avg Exp', new_col='Tier')
    create_tier(df_wt, new_col='Tier')
    # concat string for final tier name
    df_mt['Tier'] = df_mt['Tier']+df_wt['Tier']
    # df_mt['Tier'] = df_mt['Tier'].astype(int)
    return df_mt