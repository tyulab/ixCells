#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import dose_response
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
    # TODO: refactor
    # pos = "Ionis1375651"
    # neg = "Ionis676630"
    # separate _3 from _10
    df = df[['Experiment Name', 'Position', 'SampleName', 'ASO Microsynth ID', 'Exp', 'Test plate #']]# .dropna(subset=['Exp'])
    controls = df['SampleName'].str.contains('Ionis|Naive')
    df[controls].to_csv("output/controls.csv", index=False)
    control_10 = df[controls & df['SampleName'].str.contains('_10')].reset_index(drop=True)
    control_3 = df[controls & df['SampleName'].str.contains('_3')].reset_index(drop=True)
    # print(control_10)
    control_10.to_csv("output/controls_10.csv", index=False)
    control_3.to_csv("output/controls_3.csv", index=False)
    hist(control_10, "hist_10", title='10mmol Plates')
    hist(control_3, "hist_3", title='3mmol Plates')

# TODO: parameter for title of plot
def hist(control, file_prefix="hist", title='Plates'):
    # export df used for histogram
    control['SampleName'] = control['SampleName'].str.split('_').str[0]
    group_samples = control.groupby('SampleName')
    pos_df = group_samples.get_group('Ionis1375651')[['Exp', 'Test plate #']]
    try:
        neg_df = group_samples.get_group('Ionis676630')[['Exp', 'Test plate #']]
    except:
        neg_df = None

    step = 9
    layout = (3,3)
    plates = control['Test plate #'].nunique()
    for i in range(1,plates,step):
        fig, axes = plt.subplots(sharex=True, sharey=True)
        plt.margins(x=0, y=0)
        rem = plates-i+1
        axes = pos_df[pos_df['Test plate #'].between(i, i+step-1)].hist( \
            column='Exp',by=pos_df['Test plate #'], layout=layout, \
            sharex=True, sharey=True, figsize=(6,7), xrot=0, alpha=0.5, label='Ionis1375651', bins=5)
        axes = axes.ravel()[:min(rem,step)]

        # plt.yticks(np.linspace(0, 10, num=6, endpoint=True))
        plt.ylim(bottom=0,top=10)

        if neg_df is not None:
            plt.xlim(left=0.0, right=1.75)
            plt.xticks(np.linspace(0,2,num=5,endpoint=True))
            neg_df[neg_df['Test plate #'].between(i, i+step-1)].hist( \
                column='Exp', by=neg_df['Test plate #'], ax=axes, \
                alpha=0.5, xrot=0, label='Ionis676630',color='r', bins=5)
            plt.figlegend(['Ionis1375651', 'Ionis676630'], loc='lower right')
        else:
            plt.xlim(left=0.0, right=1.5)
            plt.xticks(np.linspace(0,1.5,num=4,endpoint=True))
            plt.figlegend(['Ionis1375651'], loc='lower right')
        plt.suptitle(title)
        # fig.tight_layout()


        for ax in axes.flatten():
            ax.xaxis.set_tick_params(labelbottom=True)
            ax.yaxis.set_tick_params(labelbottom=True)

        plt.savefig('plots/'+file_prefix+' '+str(i)+'-'+str(i+min(rem,step)-1)+'.png')
    # plt.show()

def main():
    # default path for table to read from
    table_file = "output/output.csv"
    plate_output = "plates_table.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # TODO: flag any controls across plates > 1.55
    outliers = dose_response.flag_outliers(df, col='Exp', greater_than=1.55)
    df = df[~outliers].reset_index(drop=True)
    # histogram
    control_hist(df)

    # # drop columns/nan rows and export
    # plates = drop_nan(df)
    # plates.to_csv("output/" + plate_output, index=False)



if __name__=="__main__":
    main()
