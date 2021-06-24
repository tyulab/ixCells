#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Dose response step- flag ASO if avg 10mmol exp value > avg 3mmol exp value

def flag_exp(df):
    # search for _3 or _10
    samples = df[df['SampleName'].str.contains("_10|_3")]
    # avg exp over samples
    mean = samples['Exp'].groupby(samples['SampleName']).transform('mean')
    df.at[df['SampleName'].str.contains("_10|_3"),'Avg_exp'] = mean
    # drop dupes
    # df.to_csv("output/avg_exp.csv", index=False)
    df = df.drop(columns='Exp')
    df = df.drop_duplicates(['SampleName', 'Avg_exp']).sort_values(['SampleName']).dropna(subset=['Avg_exp']).reset_index(drop=True)
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
    # new_df.insert(len(new_df.columns)-2,'Avg_exp_10',df['Avg_exp_10'])
    # print(df)
    # for i in range(0,len(df),2):
    #     # check contains _10 and _3
    #     if "_10" in df.loc[i,'SampleName'] and "_3" in df.loc[i+1,'SampleName']:
    # flag the ASO if the avg exp of _10 is > _3
    new_df['10mmol > 3mmol'] = new_df['Avg_exp_10'].gt(new_df['Avg_exp_3'])
    print(new_df)
    print(new_df['10mmol > 3mmol'].value_counts())
    new_df.to_csv("output/dose_response.csv", index=False)

def main():
    # default path for table to read from
    table_file = "output/output.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    df = df[['Experiment Name','SampleName','ASO Microsynth ID','Exp','Test plate #']]
    pd.set_option('display.max_columns', None)
    flag_exp(df)

if __name__=="__main__":
    main()