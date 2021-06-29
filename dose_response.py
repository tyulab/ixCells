#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Dose response step- flag ASO if avg 10mmol exp value > avg 3mmol exp value

def avg_exp(df):
    # average by replicate
    df["Avg Exp"] = np.nan
    # print(df['Exp_zscore'].groupby(df.index // 3).mean())
    mean = df['Exp'].groupby(df.index // 3).transform('mean')
    df.at[::3, 'Avg Exp'] = mean

    # change naive # TODO: move average of naive to new column
    df["Avg Exp Naive"] = np.nan
    # transform mean on naive groups
    df.at[df['SampleName'].str.contains("Naive", case=False), "Avg Exp Naive"] = \
        df[df['SampleName'].str.contains("Naive", case=False)].groupby(['Experiment Name','SampleName'])['Exp'].transform('mean')
    # set duplicate values to nan
    duplicate = df.duplicated(['SampleName', "Avg Exp Naive"], keep='first')
    df.loc[df['SampleName'].str.contains("Naive", case=False) & duplicate, "Avg Exp Naive"] = np.nan

    # get range(?) of each group (every 6)
    range = lambda x: abs(x.max() - x.min())
    # apply to exp column based on sample name
    # print(df.groupby(['Experiment Name','SampleName']).groups.keys())
    df['Exp_range'] = df.groupby(['Experiment Name','SampleName'])['Avg Exp'].transform(range)
    # set duplicate values to nan
    duplicate = df.duplicated(['SampleName', "Exp_range"], keep='first')
    df.loc[duplicate, 'Exp_range'] = np.nan

    print(df)

    df.to_csv("output/avg_exp.csv", index=False)


def flag_exp(df):
    # search for _3 or _10
    samples = df[df['SampleName'].str.contains("_10|_3")]
    # avg exp over samples
    mean = samples['Exp'].groupby(samples['SampleName']).transform('mean')
    df.at[df['SampleName'].str.contains("_10|_3"),'Avg_exp'] = mean
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
    print(new_df['10mmol > 3mmol'].value_counts())
    new_df.to_csv("output/dose_response.csv", index=False)

def main():
    # default path for table to read from
    table_file = "output/output.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    df = df[['Experiment Name','SampleName','ASO Microsynth ID','Exp','Test plate #']]
    pd.set_option('display.max_columns', None)

    # functions
    avg_exp(df)
    flag_exp(df)

if __name__=="__main__":
    main()