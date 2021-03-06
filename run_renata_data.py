#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# create histogram for controls with matplotlib
from utils import functions
from utils import stats
from utils.plots import r2_plot


# def control_hist(df):
#     # pos = "Ionis 1375651"
#     # neg = "Ionis 676630"
#     # separate _3 from _10
#     df = df[['Experiment Name','SampleName', 'ASO Microsynth ID', 'Avg Exp']]# .dropna(subset=['Avg Exp'])
#     controls = df['SampleName'].str.contains('Ionis', case=False)
#     df[controls].to_csv("Renata/output/renata_controls.csv", index=False)
#     control_10 = df[controls & df['SampleName'].str.contains('_10')].reset_index(drop=True)
#     control_3 = df[controls & df['SampleName'].str.contains('_3')].reset_index(drop=True)

#     control_10.to_csv("Renata/output/renata_controls_10.csv", index=False)
#     control_3.to_csv("Renata/output/renata_controls_3.csv", index=False)
#     hist(control_10, "renata_hist_10", title='Avg Exp 10mmol')
#     hist(control_3, "renata_hist_3", title='Avg Exp 3mmol')

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

    plt.savefig('Renata/plots/'+file_prefix+'.png')
    # plt.show()

# same as in dose response file
def avg_exp_zscore(df):
    # df = functions.hide_naive_control(df)
    df = df.replace(to_replace='_1$|_2$', value='', regex=True)
    mean = df.groupby(['Experiment Name', 'SampleName'])['Avg Exp'].transform('mean')
    std = df.groupby(['Experiment Name', 'SampleName'])['Avg Exp'].transform('std')
    df['Exp_std'] = std
    stats.zscore(df, 'Avg Exp', mean, std)
    return df


def std_hist(df, col='Exp std', color='b'):
    # separate into WT, MT, total
    # for i in range(len(types)):
    #     type_df = df_no_na[df_no_na['Experiment Name'].str.contains(types[i], case=False)]
    #     type_df.hist(column=col, bins=20)
    #     plt.title("Exp zscore ranges for "+types[i])
    #     plt.xlabel('Exp zscore range')
    #     # save the dfs in case
    #     type_df.to_csv("Renata/output/Exp zscore "+types[i]+".csv", index=False)
    #     plt.savefig('Renata/plots/' + 'Exp zscore hist ' + types[i] + '.png')


    # fig, axes = plt.subplots(sharex=True, sharey=True)
    plt.margins(x=0, y=0)
    df.hist(column=col, sharex=True, sharey=True, xrot=0, alpha=0.5, label='df', color=color, bins=40)
    # drop outliers?
    max = df[col].max()
    min = df[col].min()
    plt.xlim(left=min, right=max)

    # axes = axes.ravel()[:min(rem, step)]

    # plt.xticks(np.linspace(0, 10, num=6, endpoint=True))
    # plt.ylim(bottom=0, top=10)
    # plt.xticks(np.linspace(0, 2, num=5, endpoint=True))
    # plt.figlegend(['df1', 'df2'], loc='lower right')

    plt.show()

def make_std_hist():
    data_file_1 = "Renata/output/renata data.csv"
    data_file_2 = "Renata/output/output.csv"
    # store in df
    df1 = pd.read_csv(data_file_1, encoding='latin-1')
    df2 = pd.read_csv(data_file_2, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    df1_outliers = functions.flag_outliers(df1, col='Exp_std',greater_than=2.0)
    df2_outliers = functions.flag_outliers(df2, col='Exp_std', greater_than=2.0)
    df1 = df1[~df1_outliers]
    df2 = df2[~df2_outliers]

    std_hist(df1, col='Exp_std', color='b')
    std_hist(df2, col='Exp_std', color='r')

# change renata's data to fit ixcell scripts
def create_file():
    # default path for table to read from
    data_file = "Renata/data/renata data.csv"
    # store in dfrenata data.csv
    df = pd.read_csv(data_file, encoding='latin-1')
    df.columns = df.columns.str.strip()
    pd.set_option('display.max_columns', None)

    # sort
    df = df.sort_values(['ASO Microsynth ID'], ignore_index=True)
    df = avg_exp_zscore(df)

    #control_hist(df)
    df.to_csv("Renata/output/renata data.csv")

# helper function to separate values into belonging to biological replicate 1 or 2 and return as tuple of arrays
def get_replicates(df, type):
    # FIXME: move replicate info to new column
    b1 = df[df['Experiment Name'].str.contains('_1')].replace(to_replace='_1$|_2$', value='', regex=True)
    b2 = df[df['Experiment Name'].str.contains('_2')].replace(to_replace='_1$|_2$', value='', regex=True)
    s1 = pd.merge(b1, b2, how='inner', on=['Experiment Name', 'SampleName'],suffixes=('_1','_2'))
    s1.to_csv('Renata/output/renata_' + type + '.csv')
    return s1['Avg Exp_1'].to_numpy(),s1['Avg Exp_2'].to_numpy()

output_file = open("all_r_squared.tsv", "a") # tsv file of r2 values
def r_squared():
    file = "Renata/data/renata data.csv"
    # store in df
    df = pd.read_csv(file, encoding='latin-1')
    pd.set_option('display.max_columns', None)
    df = df[['Experiment Name','SampleName','Avg Exp']]
    df = functions.hide_naive_control(df)

    # keep only 10mmol
    #df = df[df['SampleName'].str.contains('_10')]
    #kq883_10 outlier
    df = df[~df['SampleName'].str.contains('KQ883_10')]

    pdf = matplotlib.backends.backend_pdf.PdfPages('Renata/data/r_squared.pdf')

    for type in ['WT','MT']:
        type_df = functions.get_type(df, type=type)
        b1, b2 = get_replicates(type_df, type)
        fig = r2_plot(b1,b2, title='renata_Exp40 ' + type, output_file=output_file)
        pdf.savefig(fig)
    pdf.close()

if __name__ == "__main__":
    create_file()
    r_squared()