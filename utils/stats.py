#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy import stats
import glob
import os

# Functions for commonly used stats and all calculations

# get absolute value z score from dataframe on column
from utils.functions import hide_naive_control, get_avg_exp

# bins to check intervals on. nomenclature: "10" means 0.1 <= x < 0.2, "90" means 0.9 <= x < inf, etc...
BINS = np.linspace(0.0,1.0,11) # right = False
names = list(map(lambda x: '0' + str(x),np.arange(0,100,10)))
names[0] = '000'
names.append('100')
# the name assigned to each bin by corresponding indices
NAMES = names
# key starting from 1, names at NAMES
DICT = dict(enumerate(NAMES, 1))



def abs_zscore(df, col='Exp'):
    # z score lambda function
    zscore = lambda x: abs((x - x.mean()) / x.std())
    # apply to exp column based on sample name]
    # TODO: remove std (for testing)
    df[col +'_std'] = df[col].groupby(df.index // 3).transform('std')
    df[col+'_zscore'] = df.groupby(['Experiment Name','SampleName'])[col].transform(zscore)

# redo ddCt normalized on neg control
# TODO: make it to apply after merging sheets, low priority
def neg_ddct(df):
    # get all values Ionis676.., group by experiment name/sample name
    ionis_neg = df[df['SampleName'].str.contains('Ionis676630|Ionis 676630', na=False)]
    # just more convenient to keep all columns
    df['Avg dCt Ionis676630'] = ionis_neg['dCt'].mean()
    # do difference
    df['ddCt'] = df['dCt'] - df['Avg dCt Ionis676630']

# redo exp based on ddct column
def calc_exp(df, col='ddCt', output_col='Exp'):
    expression = lambda x: 2 ** (-x)
    df[output_col] = df[col].transform(expression)

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

# avg zscore by sample
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

# calculate average zscore by replicate and range
def exp_zscore_range(df):
    # average by replicate
    df["Avg Exp"] = np.nan
    # print(df['Exp_zscore'].groupby(df.index // 3).mean())
    mean = df['Exp'].groupby(df.index // 3).transform('mean')
    df.loc[::3, 'Avg Exp'] = mean

    df = hide_naive_control(df)
    mean = df.groupby(['Experiment Name','SampleName'])['Exp'].transform('mean')
    std = df.groupby(['Experiment Name','SampleName'])['Exp'].transform('std')
    zscore(df, 'Avg Exp', mean, std)

    # get range(?) of each group (every 6)
    range = lambda x: abs(x.max() - x.min())
    # apply to exp column based on sample name
    # print(df.groupby(['Experiment Name','SampleName']).groups.keys())
    df['Avg Exp_zscore range'] = df.groupby(['Experiment Name', 'SampleName'])['Avg Exp_zscore'].transform(range)
    # # set duplicate values to nan
    # duplicate = df.duplicated(['SampleName', 'Avg Exp_zscore range'], keep='first')
    # df.loc[duplicate, 'Avg Exp_zscore range'] = np.nan

    # two_sample_ttest(df, 'Exp') # confirm results
    return df

# helper fn: pass MT and WT, assign to tiers
def create_tier(df, col='Avg Exp', new_col='Tier'):
    df[new_col] = np.vectorize(DICT.get)(np.digitize(df[col], BINS, right=False))


def tierlist(df):
    # get Avg Exp, group by sample
    df['Avg Exp'] = df.groupby(['Experiment Name','SampleName'])['Avg Exp'].transform('mean')
    # df.to_csv("output/avg_exp.csv", index=False)
    df_mt = df[df['Experiment Name'].str.contains('MT', case=False)].reset_index(drop=True)[['SampleName','ASO Microsynth ID','Test plate #','Avg Exp']]
    df_wt = df[df['Experiment Name'].str.contains('WT', case=False)].reset_index(drop=True)[['SampleName','ASO Microsynth ID','Test plate #','Avg Exp']]
    df = df[['SampleName','ASO Microsynth ID','Avg Exp']]
    # df['Tier'] = np.nan
    df_mt = df_mt.groupby(df_mt['SampleName']).first()
    df_mt = df_mt.rename(columns={'Avg Exp': 'MT Avg Exp'})
    df_wt = df_wt.groupby(df_wt['SampleName']).first()
    df_mt.to_csv("output/df_mt.csv")
    df_wt.to_csv("output/df_wt.csv")
    # TODO: make sure count is even
    df_mt['WT Avg Exp'] = df_wt['Avg Exp']
    create_tier(df_mt, col='MT Avg Exp', new_col='Tier')
    create_tier(df_wt, new_col='Tier')
    df_mt['Tier'] = df_mt['Tier']+df_wt['Tier']
    # df_mt['Tier'] = df_mt['Tier'].astype(int)
    return df_mt
    # return df.sort_values(['Tier', 'ASO Microsynth ID'], na_position='last')


# old helper fn to pass MT and WT, assign to tiers
def rank_tier(mt, wt):
    tier = np.nan
    if mt < 0.2:  # tier 1
        if wt > 0.8:
            tier = '1A'
        elif wt > 0.7:
            tier = '1B'
        elif wt > 0.6:
            tier = '1C'
    elif mt < 0.3:  # tier 2
        if wt > 0.6:
            tier = '2A'
        elif wt > 0.5:
            tier = '2B'
        elif wt > 0.4:
            tier = '2C'
        elif wt > 0.3:
            tier = '2D'
    elif mt < 0.4:  # tier 3
        if wt > 0.7:
            tier = '3A'
        elif wt > 0.6:
            tier = '3B'
    elif mt < 0.5 and wt > 0.7:
        tier = '3C'
    return tier


def tierlist2(df):
    # get Avg Exp, group by sample
    df['Tier'] = np.nan
    df['Avg Exp'] = df.groupby(['Experiment Name','SampleName'])['Avg Exp'].transform('mean')
    df.to_csv("output/avg_exp.csv")
    df_mt = df[df['Experiment Name'].str.contains('MT', case=False)].reset_index(drop=True)
    df_wt = df[df['Experiment Name'].str.contains('WT', case=False)].reset_index(drop=True)
    df_mt = df_mt['Avg Exp'].groupby(df_mt['SampleName']).mean()
    df_wt = df_wt['Avg Exp'].groupby(df_wt['SampleName']).mean()
    df_mt.to_csv("output/df_mt.csv")
    df_wt.to_csv("output/df_wt.csv")
    samples = df['SampleName'].unique()
    count = 0
    for sample in samples:
        mt_avg_exp = get_avg_exp(df_mt, sample)
        wt_avg_exp = get_avg_exp(df_wt, sample)
        if mt_avg_exp is not None and wt_avg_exp is not None:
            count += 1
            tier = rank_tier(mt_avg_exp, wt_avg_exp)
            df.loc[df['SampleName'].str.contains(sample),'Tier'] = tier
    print('count = ' + str(count))
    return df.sort_values(['Tier', 'ASO Microsynth ID'], na_position='last')

# compare 10mmol vs 3mmol
def flag_exp(df):
    # search for _3 or _10
    samples = df[df['SampleName'].str.contains("_10|_3")]
    # avg exp over samples
    mean = samples['Exp'].groupby(samples['SampleName']).transform('mean')
    df.loc[df['SampleName'].str.contains("_10|_3"),'Avg_exp'] = mean
    # drop dupes and exp column
    df = df.drop(columns='Exp')
    df = df.drop_duplicates(['SampleName', 'Avg_exp']).sort_values(['SampleName']).dropna(subset=['Avg_exp']).reset_index(drop=True)
    # separate _10 and _3 and merge
    _10 = df[df['SampleName'].str.contains("_10")].copy()
    _10['SampleName'] = _10['SampleName'].str.replace("_10","")
    _3 = df[df['SampleName'].str.contains("_3")].copy()
    _3['SampleName'] = _3['SampleName'].str.replace("_3", "")
    new_df = _10.merge(_3, on=['Experiment Name', 'SampleName', 'Test plate #'], suffixes=("_10","_3"))
    cols = new_df.columns.tolist()
    cols[2:] = sorted(cols[2:])
    new_df = new_df[cols]
    plates = new_df['Test plate #']
    new_df.drop(columns=['Test plate #'], inplace=True)
    new_df.insert(2,'Test plate #', plates)
    # make flag col based on two columns
    new_df['10mmol > 3mmol'] = new_df['Avg_exp_10'].gt(new_df['Avg_exp_3'])
    # print(new_df['10mmol > 3mmol'].value_counts())
    new_df.to_csv("output/dose_response.csv", index=False)