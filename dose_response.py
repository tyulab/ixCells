#!/usr/bin/env python
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import glob
import os

# Dose response step- flag ASO if avg 10mmol exp value > avg 3mmol exp value

# Calculate mean exp for each biological replicate
# TODO: zscore
# TODO: take out naive/controls
def avg_exp(df):
    # average by replicate
    df["Avg Exp"] = np.nan
    # print(df['Exp_zscore'].groupby(df.index // 3).mean())
    mean = df['Exp'].groupby(df.index // 3).transform('mean')
    df.loc[::3, 'Avg Exp'] = mean

    # # change naive # TODO: move average of naive to new column
    # df["Avg Exp Naive"] = np.nan
    # # transform mean on naive groups
    # df.loc[df['SampleName'].str.contains("Naive", case=False), "Avg Exp Naive"] = \
    #     df[df['SampleName'].str.contains("Naive", case=False)].groupby(['Experiment Name','SampleName'])['Exp'].transform('mean')
    # # set duplicate values to nan
    # duplicate = df.duplicated(['SampleName', "Avg Exp Naive"], keep='first')
    # df.loc[df['SampleName'].str.contains("Naive", case=False) & duplicate, "Avg Exp Naive"] = np.nan


    # TODO: zscore on all 6 values?
    df = hide_naive_control(df)
    mean = df.groupby(['Experiment Name','SampleName'])['Exp'].transform('mean')
    std = df.groupby(['Experiment Name','SampleName'])['Exp'].transform('std')
    zscore(df, 'Avg Exp', mean, std)

    # get range(?) of each group (every 6)
    range = lambda x: abs(x.max() - x.min())
    # apply to exp column based on sample name
    # print(df.groupby(['Experiment Name','SampleName']).groups.keys())
    df['Avg Exp_zscore range'] = df.groupby(['Experiment Name', 'SampleName'])['Avg Exp_zscore'].transform(range)
    # set duplicate values to nan
    duplicate = df.duplicated(['SampleName', 'Avg Exp_zscore range'], keep='first')
    df.loc[duplicate, 'Avg Exp_zscore range'] = np.nan

    two_sample_ttest(df, 'Exp') # confirm results

    df.to_csv("output/avg_exp_renormalized_to_neg.csv", index=False)
    return df


def zscore(df, col, mean, std):
    # zscore = lambda x: abs((x - mean) / std)
    # df[col+'_zscore'] = df.groupby(['Experiment Name','SampleName'])[col].transform(zscore)
    df[col+'_zscore'] = df[col].sub(mean).div(std)# .abs()
    # for i in range(len(df)):
    #     if pd.notna(df.loc[i, col]):
    #     df.loc[i,col+'_zscore'] = abs((df.loc[i, col] - mean.loc[i]))

def two_sample_ttest(df, col):
    df[col+' t_stat'] = df[col + ' p_value'] = df[col + ' dof'] = np.nan
    # iterate over every 6 rows
    for i in range(0, len(df), 6):
        t, p, dof = welch_ttest(df[i:i+3][col], df[i+3:i+6][col])
        df.loc[i,[col+' t_stat', col + ' p_value', col + ' dof']] = [t, p, dof]
    # print(df)


def welch_ttest(x, y):
    ## Welch-Satterthwaite Degrees of Freedom ##
    # print(x)
    # print(y)
    t, p = stats.ttest_ind(x, y, equal_var=False, nan_policy='omit')
    return t, p, welch_dof(x, y)

    # print("\n",
    #       f"Welch's t-test= {t:.4f}", "\n",
    #       f"p-value = {p:.4f}", "\n",
    #       f"Welch-Satterthwaite Degrees of Freedom= {dof:.4f}")

def welch_dof(x, y):
    dof = (x.var() / x.size + y.var() / y.size) ** 2 / ((x.var() / x.size) ** 2 / (x.size - 1) + (y.var() / y.size) ** 2 / (y.size - 1))
    return dof


def hide_naive_control(df):
    df = df[~df['SampleName'].str.contains('Naive|Ionis1375651|Ionis676630', case=False)].copy().reset_index(drop=True)
    return df

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

# set outliers above a certain z score threshold to na (preserve shape)
def clear_outliers(df, zscore_col='Exp_zscore', exp_col='Exp', greater_than=1.8):
    count = df[df[zscore_col] > greater_than].count()[0]
    df.loc[df[zscore_col] > greater_than, exp_col] = np.nan
    print("%s: %d greater than %s" % (zscore_col, count, greater_than))
    return df

# helper fn. count outliers and return
def flag_outliers(df, zscore_col='Exp_zscore', greater_than=2):
    count = df[df[zscore_col] > greater_than].count()[0]
    # df[zscore_col+'<'+greater_than] = True
    # df.loc[df[zscore_col] > greater_than, exp_col] = False
    print("%s: %d greater than %s" % (zscore_col, count, greater_than))
    return df[zscore_col] > greater_than

def type_hist(df, col):
    # separate into WT, MT, total
    df_no_na = df.dropna(subset=[col]).reset_index(drop=True)
    types = ['WT','MT','total']
    for i in range(len(types)):
        type_df = df_no_na[df_no_na['Experiment Name'].str.contains(types[i], case=False)]
        type_df.hist(column=col, bins=20)
        plt.title("Exp zscore ranges for "+types[i])
        plt.xlabel('Exp zscore range')
        # save the dfs in case
        type_df.to_csv("output/Exp zscore "+types[i]+".csv", index=False)
        plt.savefig('plots/' + 'Exp zscore hist ' + types[i] + '.png')
    # plt.show()

def rank_tier(mt, wt):
    tier = np.nan
    if mt < 0.2:  # tier 1
        if wt > 0.7:
            tier = '1A'
        elif wt > 0.6:
            tier = '1B'
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

def tierlist(df):
    df['Tier'] = np.nan
    df_mt = df[df['Experiment Name'].str.contains('MT', case=False)].reset_index(drop=True)
    df_wt = df[df['Experiment Name'].str.contains('WT', case=False)].reset_index(drop=True)
    # print(df_mt)
    for row in range(len(df_mt)):
        sample = df_mt.loc[row,'SampleName']
        mt_avg_exp = df_mt.loc[row,'Avg Exp']
        wt_avg_exp = df_wt.loc[row,'Avg Exp']
        tier = rank_tier(mt_avg_exp, wt_avg_exp)
        df.loc[df['SampleName'].str.contains(sample),'Tier'] = tier
    return df.sort_values(['Tier', 'SampleName'], na_position='last')

def main():
    # default path for table to read from
    table_file = "output/output.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    df = df[['Experiment Name','SampleName','ASO Microsynth ID','Exp_zscore','Exp','Test plate #']]
    pd.set_option('display.max_columns', None)

    # functions
    # df = clear_outliers(df)
    df = avg_exp(df)
    # flag_exp(df)

    type_hist(df, 'Avg Exp_zscore range')

    # flag >1.55 zscore range
    df = df.dropna(subset=['Avg Exp_zscore range']).reset_index(drop=True)
    ranges = flag_outliers(df, zscore_col='Avg Exp_zscore range', greater_than=1.55)
    df = df[~ranges].reset_index(drop=True)
    df = tierlist(df)
    df.to_csv("output/tiers_ungrouped.csv")
    grouped_samples = df.groupby('SampleName').first().sort_values(['Tier', 'SampleName'], na_position='last')
    grouped_samples = grouped_samples[['ASO Microsynth ID','Test plate #','Avg Exp','Tier']]
    grouped_samples.to_csv("output/tiers.csv")
    # drop ranges



if __name__ == "__main__":
    main()