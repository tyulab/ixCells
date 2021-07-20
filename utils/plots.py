#!/usr/bin/env python
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import colors
from tqdm import tqdm
import glob
import os
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

def type_hist(df, col):
    # separate into WT, MT, total
    df_no_na = df.dropna(subset=[col]).reset_index(drop=True)
    types = ['WT','MT','total']
    for i in range(len(types)):
        type_df = df_no_na[df_no_na['Experiment Name'].str.contains(types[i], case=False)]
        type_df.hist(column=col, bins=20)
        plt.title("Exp zscore ranges for "+types[i])
        plt.xlabel('Exp zscore range')
        # save the dfs in case
        type_df.to_csv("output/Exp zscore "+types[i]+".csv", index=False)
        plt.savefig('plots/' + 'Exp zscore hist ' + types[i] + '.png')
    # plt.show()

# create histogram for controls with matplotlib
def control_hist(df):
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
    if not control_10.empty:
        hist_plates(control_10, "hist_10", title='10mmol Plates')
    if not control_3.empty:
        hist_plates(control_3, "hist_3", title='3mmol Plates')

# create histograms by plates
def hist_plates(control, file_prefix="hist", title='Plates'):
    # export df used for histogram
    control['SampleName'] = control['SampleName'].str.split('_').str[0]
    group_samples = control.groupby('SampleName')
    # print(group_samples.groups)
    pos_df = group_samples.get_group('Ionis1375651')
    pos_df = pos_df[['Exp', 'Test plate #']]
    try:
        neg_df = group_samples.get_group('Ionis676630')[['Exp', 'Test plate #']]
    except:
        neg_df = None

    step = 9
    layout = (3,3)
    plates = control['Test plate #'].nunique()
    for i in tqdm(range(1,plates+1,step)):
        # fig, axes = plt.subplots(sharex=True, sharey=True)
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

# helper function to create scatter plot for the tiers df
def tier_scatter(df, title='Tiers by MT/WT'):
    # Create nums for MT and WT
    df['MT Avg Exp'] = df['MT Avg Exp']*100
    df['WT Avg Exp'] = df['WT Avg Exp']*100

    cmap = colors.ListedColormap(['r', 'g', 'b', 'c'])
    bounds = [0, 10, 20]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    ax1 = df.plot.scatter(x='WT Avg Exp', y = 'MT Avg Exp') # alpha=0.5

    plt.title(title)
    plt.savefig('plots/scatter_uncropped.pdf')

    df = df[df['MT Avg Exp'] <= 110]
    df = df[df['WT Avg Exp'] <= 110]
    # ax1 = df.plot.scatter(x='WT Avg Exp', y = 'MT Avg Exp', c=df['MT Avg Exp'], cmap=cmap)
    plt.xlim(left=0.0, right=100)
    plt.ylim(bottom=0, top=100)
    plt.xticks(np.linspace(0,110,num=12,endpoint=True))
    plt.yticks(np.linspace(0,110,num=12,endpoint=True))
    plt.grid(True)

    plt.savefig('plots/scatter.pdf')
    plt.show()


# helper function to create scatter plot for the tiers df
def tier_hist(df):
    # Separate the tiers for MT and WT
    mt_tier = df['Tier'].str.slice(stop=3)
    wt_tier = df['Tier'].str.slice(start=3)
    # print(mt_tier)

    # # export df used for histogram
    # control['SampleName'] = control['SampleName'].str.split('_').str[0]
    # group_samples = control.groupby('SampleName')
    # pos_df = group_samples.get_group('Ionis1375651')[['Exp', 'Test plate #']]
    # try:
    #     neg_df = group_samples.get_group('Ionis676630')[['Exp', 'Test plate #']]
    # except:
    #     neg_df = None
    #
    # step = 9
    # layout = (3, 3)
    # plates = control['Test plate #'].nunique()
    # for i in range(1, plates + 1, step):
    #     # fig, axes = plt.subplots(sharex=True, sharey=True)
    #     plt.margins(x=0, y=0)
    #     rem = plates - i + 1
    #     axes = pos_df[pos_df['Test plate #'].between(i, i + step - 1)].hist_plates( \
    #         column='Exp', by=pos_df['Test plate #'], layout=layout, \
    #         sharex=True, sharey=True, figsize=(6, 7), xrot=0, alpha=0.5, label='Ionis1375651', bins=5)
    #     axes = axes.ravel()[:min(rem, step)]
    #
    #     # plt.yticks(np.linspace(0, 10, num=6, endpoint=True))
    #     plt.ylim(bottom=0, top=10)
    #
    #     if neg_df is not None:
    #         plt.xlim(left=0.0, right=1.75)
    #         plt.xticks(np.linspace(0, 2, num=5, endpoint=True))
    #         neg_df[neg_df['Test plate #'].between(i, i + step - 1)].hist_plates( \
    #             column='Exp', by=neg_df['Test plate #'], ax=axes, \
    #             alpha=0.5, xrot=0, label='Ionis676630', color='r', bins=5)
    #         plt.figlegend(['Ionis1375651', 'Ionis676630'], loc='lower right')
    #     else:
    #         plt.xlim(left=0.0, right=1.5)
    #         plt.xticks(np.linspace(0, 1.5, num=4, endpoint=True))
    #         plt.figlegend(['Ionis1375651'], loc='lower right')
    #     plt.suptitle(title)
    #     # fig.tight_layout()
    #
    #     for ax in axes.flatten():
    #         ax.xaxis.set_tick_params(labelbottom=True)
    #         ax.yaxis.set_tick_params(labelbottom=True)
    #
    #     plt.savefig('plots/' + file_prefix + ' ' + str(i) + '-' + str(i + min(rem, step) - 1) + '.png')
    # # plt.show()


def r2_plot(y_test,y_predicted, title):
    fig, ax = plt.subplots()
    ax.scatter(y_test, y_predicted)
    # ax.plot([y_test.min(), y_test.max()], [y_predicted.min(), y_predicted.max()], 'k--', lw=4)
    ax.set_xlabel('B1')
    ax.set_ylabel('B2')
    # regression line
    y_test, y_predicted = y_test.reshape(-1, 1), y_predicted.reshape(-1, 1)
    ax.plot(y_test, LinearRegression().fit(y_test, y_predicted).predict(y_test), color='red')
    ax.set_title(title+' R2: ' + str(r2_score(y_test, y_predicted)))
    plt.savefig('plots/'+title+'_R2.pdf')
    plt.show()