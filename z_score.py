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

    # change naive # TODO: fix str.contains
    # transform mean on naive groups
    df.at[df['SampleName'].str.contains("Naive", case=False),"Avg_zscore"] = df[df['SampleName'].str.contains("Naive", case=False)].groupby('SampleName')['Exp_zscore'].transform('mean')
    # set duplicate values to nan
    duplicate = df.duplicated(['SampleName', 'ASO Microsynth ID', 'Avg_zscore'], keep='first')
    df.loc[df['SampleName'].str.contains("Naive", case=False) & duplicate, 'Avg_zscore'] = np.nan

# avg zscore by plate
def avg_plates(df):
    # get average by plate value
    df["Avg_plate"] = np.nan
    mean = df['Exp_zscore'].groupby(df['Test plate #']).transform('mean')
    df["Avg_plate"] = mean
    # set duplicate values to nan
    duplicate = df.duplicated(['Test plate #', 'Avg_plate'],keep='first')
    df.at[duplicate, 'Avg_plate'] = np.nan

# get list of all csvs
def get_csv(path, plate_file):
    # get all csvs in path, ignore summary data file, plate info, and output folder
    csv = [f for f in glob.iglob(os.path.join(path, "**/*.csv"),recursive=True) \
           if not 'summary data' in os.path.basename(f) and not plate_file in os.path.basename(f)
           and not 'output' in os.path.dirname(f)]
    return csv

# read plates and add column to output
def assign_plates(df, file="ixCells_Round 1_2021-06-22_TN09_551ASOs_plate id adjusted.csv"):
    plates = pd.read_csv(file, encoding='latin-1')
    plates.columns = plates.columns.str.strip()
    plates = plates[['Test plate #', 'ASO Microsynth ID']]
    df["Test plate #"] = np.nan
    for i in range(len(df)):
        aso = str(df.loc[i]['ASO Microsynth ID']).split("_")[0]
        if aso in plates['ASO Microsynth ID'].values:
            df.at[i,'Test plate #'] = plates[plates['ASO Microsynth ID'] == aso]['Test plate #']

def main():
    # TODO: parser?
    # specify path to folder containing all csvs and plate sheet + output to be ignored
    path = os.getcwd()
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
            # df['SampleName'] = df['SampleName'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')
            # df['ASO Microsynth ID'] = df['ASO Microsynth ID'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')

            df_updated = df.replace(to_replace='^Na.*ve', value='Naive', regex=True)
            df_list.append(df_updated)
    concat_df = pd.concat(df_list, ignore_index=True)
    concat_df = concat_df[['Experiment Name', 'SampleName','ASO Microsynth ID', 'ddCt','Exp']]
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