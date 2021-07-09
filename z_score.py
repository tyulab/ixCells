#!/usr/bin/env python
import pandas as pd
import numpy as np
import glob
import os
import re
# Loads csvs, gets zscore, average zscore by sample and plate, and puts to output csv

# set outliers above a certain z score threshold to na (preserve shape)
def clear_outliers(df, zscore_col='Exp_zscore', exp_col='Exp', greater_than=1.8):
    count = df[df[zscore_col] > greater_than].count()[0]
    df.loc[df[zscore_col] > greater_than, exp_col] = np.nan
    print("%s: %d greater than %s" % (zscore_col, count, greater_than))
    return df

# helper fn. count outliers and return
def flag_outliers(df, col='Exp_zscore', greater_than=2):
    count = df[df[col] > greater_than].count()[0]
    # df[zscore_col+'<'+greater_than] = True
    # df.loc[df[col] > greater_than, exp_col] = False
    print("%s: %d greater than %s" % (col, count, greater_than))
    return df[col] > greater_than

# get absolute value z score from dataframe on column
def abs_zscore(df, col='Exp'):
    # z score lambda function
    zscore = lambda x: abs((x - x.mean()) / x.std())
    # apply to exp column based on sample name]
    # TODO: remove std (for testing)
    df[col +'_std'] = df[col].groupby(df.index // 3).transform('std')
    df[col+'_zscore'] = df.groupby(['Experiment Name','SampleName'])[col].transform(zscore)

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

# avg zscore by plate
def avg_plates(df):
    # fixed: includes naive and control into calculations
    # sort on test plate and sample
    df["Test plate #"] = pd.to_numeric(df["Test plate #"], errors='coerce')
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

# redo ddCt normalized on neg control
# TODO: make it to apply after merging sheets, but it's a hassle
def neg_ddct(df):
    # get all values Ionis676.., group by experiment name/sample name
    ionis_neg = df[df['SampleName'].str.contains('Ionis676630')]
    # just more convenient to keep all columns
    df['Avg dCt Ionis676630'] = ionis_neg['dCt'].mean()
    # do difference
    df['ddCt'] = df['dCt'] - df['Avg dCt Ionis676630']

# redo exp based on ddct column
def redo_exp(df, col='ddCt', output_col='Exp'):
    expression = lambda x: 2 ** (-x)
    df[output_col] = df[col].transform(expression)

# read plates and add column to output
def assign_plates(df, file="ixCells_Round 1_2021-06-22_TN09_551ASOs_plate id adjusted.csv"):
    plates = pd.read_csv(file, encoding='latin-1')
    plates.columns = plates.columns.str.strip()
    plates = plates[['Test plate #', 'ASO Microsynth ID']]
    df["Test plate #"] = np.nan
    for i in range(len(df)): # try to find _P[0-9]+
        # sample = str(df.loc[i]['SampleName']).split("_")
        find_plate_no = re.search('_P\d+', str(df.loc[i]['SampleName'])) # search for plate number
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

# get list of all csvs
def get_csv(path, plate_file):
    # get all csvs in path holding data, ignore files that don't start with ASO
    csv = [f for f in glob.iglob(os.path.join(path, "**/*.csv"),recursive=True) \
           if os.path.basename(f).startswith('ASO')]
    return csv

def main():
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
            neg_ddct(df)
            df = df.replace(to_replace='^Na.*ve', value='Naive', regex=True)
            df_list.append(df)

    df = pd.concat(df_list, ignore_index=True)
    # pd.set_option('display.max_rows', df.shape[0]+1)
    # calculate ddCt
    # Calculate exp
    redo_exp(df)
    redo_exp(df, col='ddCt from Naive', output_col='Exp from Naive')
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
    # concat_df.sort_values(['SampleName', 'Experiment Name'], ignore_index=True, inplace=True)
    print(df)
    # export to output file
    df.to_csv("output/"+output_file, index=False)


if __name__ == "__main__":
    main()
