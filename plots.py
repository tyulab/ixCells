#!/usr/bin/env python

# make histogram for +/- control and output file with plate values using output file

# drop the unnecessary rows and send to csv to view plates only
from run import create_control_hist


def drop_nan(df):
    naive_df = df.dropna(subset=['Avg_naive'])
    naive_df = naive_df[['Test plate #','Experiment Name', 'Avg_naive']]
    print(naive_df)
    new_df = df.dropna(subset=['Avg_plate'])
    new_df = new_df[['Test plate #', 'Avg_plate', 'Min_plate', 'Max_plate']]
    new_df = new_df.merge(naive_df, on='Test plate #').drop_duplicates()
    print(new_df)
    return new_df


if __name__=="__main__":
    create_control_hist()
