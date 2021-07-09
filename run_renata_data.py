#!/usr/bin/env python
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import glob
import os


# create histogram for controls with matplotlib
def control_hist(df):
    # pos = "Ionis 1375651"
    # neg = "Ionis 676630"
    # separate _3 from _10
    df = df[['Experiment Name','SampleName', 'ASO Microsynth ID', 'Avg Exp']]# .dropna(subset=['Avg Exp'])
    controls = df['SampleName'].str.contains('Ionis|Naive')
    df[controls].to_csv("output/renata_controls.csv", index=False)
    control_10 = df[controls & df['SampleName'].str.contains('_10')].reset_index(drop=True)
    control_3 = df[controls & df['SampleName'].str.contains('_3')].reset_index(drop=True)

    control_10.to_csv("output/renata_controls_10.csv", index=False)
    control_3.to_csv("output/renata_controls_3.csv", index=False)
    hist(control_10, "renata_hist_10", title='Avg Exp 10mmol')
    hist(control_3, "renata_hist_3", title='Avg Exp 3mmol')

# TODO: parameter for title of plot
def hist(control, file_prefix="hist", title='Plates'):
    # export df used for histogram
    control['SampleName'] = control['SampleName'].str.split('_').str[0]
    group_samples = control.groupby('SampleName')
    pos_df = group_samples.get_group('Ionis 1375651')[['Avg Exp']]
    try:
        neg_df = group_samples.get_group('Ionis 676630')[['Avg Exp']]
    except:
        neg_df = None

    fig, axes = plt.subplots(sharex=True, sharey=True)
    plt.margins(x=0, y=0)
    axes = pos_df.hist( \
        column='Avg Exp', \
        sharex=True, sharey=True, figsize=(6,7), xrot=0, alpha=0.5, label='Ionis1375651', bins=5)
    axes = axes.ravel()[:]

    # plt.yticks(np.linspace(0, 10, num=6, endpoint=True))
    plt.ylim(bottom=0,top=10)

    if neg_df is not None:
        plt.xlim(left=0.0, right=1.75)
        plt.xticks(np.linspace(0,2,num=5,endpoint=True))
        neg_df.hist( \
            column='Avg Exp', ax=axes, \
            alpha=0.5, xrot=0, label='Ionis676630',color='r', bins=5)
        plt.figlegend(['Ionis1375651', 'Ionis676630'], loc='lower right')
    else:
        plt.xlim(left=0.0, right=1.5)
        plt.xticks(np.linspace(0,1.5,num=4,endpoint=True))
        plt.figlegend(['Ionis1375651'], loc='lower right')
        # fig.tight_layout()

    plt.title(title)
    for ax in axes.flatten():
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.yaxis.set_tick_params(labelbottom=True)

    plt.savefig('plots/'+file_prefix+'.png')
    # plt.show()

# change renata's data to fit ixcell scripts
def main():
    # default path for table to read from
    data_file = "data/renata data.csv"
    # store in df
    df = pd.read_csv(data_file, encoding='latin-1')
    df.columns = df.columns.str.strip()
    pd.set_option('display.max_columns', None)

    # sort
    df = df.sort_values(['ASO Microsynth ID'], ignore_index=True)
    print(df)
    df.to_csv("output/renata data.csv")

    control_hist(df)



if __name__ == "__main__":
    main()