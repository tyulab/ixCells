#!/usr/bin/env python
import pandas as pd
import numpy as np
import glob
import os
# Loads csvs, gets zscore, average zscore by sample and plate, and puts to output csv

# get absolute value z score from dataframe
def abs_zscore(df):
    # z score lambda function
    zscore = lambda x: abs((x - x.mean()) / x.std())
    # apply to exp column based on sample name
    df['Exp_zscore'] = df.groupby('SampleName')['Exp'].transform(zscore)

# avg zscore by sample
def avg_zscore(df):
    # average by replicate
    df["Avg_zscore"] = np.nan
    # print(df['Exp_zscore'].groupby(df.index // 3).mean())
    mean = df['Exp_zscore'].groupby(df.index // 3).transform('mean')
    df.at[::3,'Avg_zscore'] = mean

    # change naive # TODO: move average of naive to new column
    df["Avg_naive"] = np.nan
    # transform mean on naive groups
    df.at[df['SampleName'].str.contains("Naive", case=False),"Avg_naive"] = df[df['SampleName'].str.contains("Naive", case=False)].groupby('SampleName')['Exp_zscore'].transform('mean')
    # set duplicate values to nan
    duplicate = df.duplicated(['SampleName', 'ASO Microsynth ID', 'Avg_naive'], keep='first')
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
    df.at[duplicate, ['Avg_plate', "Min_plate", "Max_plate"]] = np.nan



# read plates and add column to output
def assign_plates(df, file="ixCells_Round 1_2021-06-22_TN09_551ASOs_plate id adjusted.csv"):
    plates = pd.read_csv(file, encoding='latin-1')
    plates.columns = plates.columns.str.strip()
    plates = plates[['Test plate #', 'ASO Microsynth ID']]
    df["Test plate #"] = np.nan
    for i in range(len(df)):
        sample = str(df.loc[i]['SampleName']).split("_")
        # add plate #s to naive and control
        if "Ionis" in sample[0] or "Naive" in sample[0]:  # check if control or naive
            # get plate number from end eg _P##
            df.at[i, 'Test plate #'] = int(sample[-1][1:])
        else:
            aso = str(df.loc[i]['ASO Microsynth ID']).split("_")[0]
            if aso in plates['ASO Microsynth ID'].values:
                df.at[i,'Test plate #'] = int(plates[plates['ASO Microsynth ID'] == aso]['Test plate #'])

# get list of all csvs
def get_csv(path, plate_file):
    # get all csvs in path, ignore summary data file, plate info, and output folder
    csv = [f for f in glob.iglob(os.path.join(path, "**/*.csv"),recursive=True) \
           if not 'summary data' in os.path.basename(f) and not plate_file in os.path.basename(f)
           and not 'output' in os.path.dirname(f)]
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
    # read all csvs into dataframes and concatenate
    for f in csv:
        if plate_file not in f:
            df = pd.read_csv(f, encoding='latin-1')
            df.columns = df.columns.str.strip()
            df_updated = df.replace(to_replace='^Na.*ve', value='Naive', regex=True)
            df_list.append(df_updated)
    concat_df = pd.concat(df_list, ignore_index=True)
    # drop irrelevant columns
    concat_df = concat_df[['Experiment Name','Position', 'SampleName','ASO Microsynth ID', 'ddCt','Exp']]
    # pd.set_option('display.max_rows', df.shape[0]+1)
    # Calculate z scores for Exp
    abs_zscore(concat_df)
    # Average z scores
    avg_zscore(concat_df)
    assign_plates(concat_df, plate_file)
    avg_plates(concat_df) # can separate this part later?
    print(concat_df)
    # export to output file
    concat_df.to_csv("output/"+output_file, index=False)

if __name__=="__main__":
    main()