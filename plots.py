#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# make histogram for +/- control and output file with plate values using output file

# drop the unnecessary rows and send to csv to view plates only
def drop_nan(df):
    naive_df = df.dropna(subset=['Avg_naive'])
    naive_df = naive_df[['Test plate #','Experiment Name', 'Avg_naive']]
    print(naive_df)
    new_df = df.dropna(subset=['Avg_plate'])
    new_df = new_df[['Test plate #', 'Avg_plate', 'Min_plate', 'Max_plate']]
    new_df = new_df.merge(naive_df, on='Test plate #').drop_duplicates()
    print(new_df)
    return new_df

# create histogram for controls with matplotlib
def control_hist(df):
    # pos = "Ionis1375651"
    # neg = "Ionis676630"
    control = df[df['SampleName'].str.contains('Ionis1375651|Ionis676630')].reset_index(drop=True)
    # export df used for histogram
    control = control[['SampleName', 'Exp', 'Test plate #']].dropna(subset=['Exp'])
    control.to_csv("output/controls.csv", index=False)
    control['SampleName'] = control['SampleName'].str.split('_').str[0]
    group_samples = control.groupby('SampleName')
    pos_df = group_samples.get_group('Ionis1375651')[['Exp', 'Test plate #']]
    neg_df = group_samples.get_group('Ionis676630')[['Exp', 'Test plate #']]
    # TODO: fix getting plates
    step = 8
    for i in range(1,control['Test plate #'].nunique(),step):
        axes = pos_df[pos_df['Test plate #'].between(i, i+step-1)].hist(column='Exp',by=pos_df['Test plate #'], layout=(2,-(-step//2)), alpha=0.5, label='Ionis1375651', bins=5)
        neg_df[neg_df['Test plate #'].between(i, i+step-1)].hist(column='Exp', by=neg_df['Test plate #'], ax=axes.ravel()[:step], alpha=0.5, label='Ionis676630',color='r', bins=5)
        plt.suptitle('Plates')
        plt.figlegend(['Ionis1375651', 'Ionis676630'], loc='upper right')
        plt.savefig('plots/hist'+str(i)+'-'+str(i+step-1)+'.png')
    # plt.show()

def main():
    # default path for table to read from
    table_file = "output/output.csv"
    plate_output = "plates_table.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # histogram
    # control_hist(df)

    # # drop columns/nan rows and export
    plates = drop_nan(df)
    plates.to_csv("output/" + plate_output, index=False)



if __name__=="__main__":
    main()
