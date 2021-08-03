#!/usr/bin/env python
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.transforms import Affine2D
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from tqdm import tqdm
import glob
import os
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

from utils.stats import *

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

    # cmap = colors.ListedColormap(['r', 'g', 'b', 'c'])
    # bounds = [0, 10, 20]
    # norm = colors.BoundaryNorm(bounds, cmap.N)

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
    ax.plot(y_test, LinearRegression().fit(y_test, y_predicted).predict(y_test), color='red', linewidth=0.75)
    ax.set_title(title+' R2: ' + str(r2_score(y_test, y_predicted)))
    plt.savefig('plots/'+title+'_R2.pdf')
    plt.show()

def box_plot(df):
    # plt.subplots(1, 1, figsize=(15, 10))
    # separate by plates
    df.loc[df['Experiment Name'].str.contains('MT'), 'Experiment Name'] = 'MT'
    df.loc[df['Experiment Name'].str.contains('Total'), 'Experiment Name'] = 'Total'
    df.loc[df['Experiment Name'].str.contains('WT'), 'Experiment Name'] = 'WT'
    plates = df['Test plate #'].nunique()
    for plate_no in range(1,plates+1):
        # fig, ax = plt.subplots()
        plate_group = df[df['Test plate #'] == plate_no]
        plate_group = plate_group.dropna(subset=['Avg Exp']).reset_index(drop=True)

        # plate_group.boxplot(column='Avg Exp', by=['SampleName','Experiment Name'], rot=45, figsize=(20, 10))

        plt.xticks(ha='right')
        plt.ylabel('Average Expression')
        plt.title('Plate ' + str(plate_no) + ' Avg Exp')
        plt.suptitle(None)

        # seaborn
        # sns.set_theme(style="whitegrid")
        ax = sns.catplot(x="SampleName", y="Avg Exp", hue="Experiment Name",
                         data=plate_group, kind="box", height=10, aspect=15/10)
        plt.show()

def error_plot(df):
    # plt.subplots(1, 1, figsize=(15, 10))
    # separate by plates
    df.loc[df['Experiment Name'].str.contains('MT'), 'Experiment Name'] = 'MT'
    df.loc[df['Experiment Name'].str.contains('Total'), 'Experiment Name'] = 'Total'
    df.loc[df['Experiment Name'].str.contains('WT'), 'Experiment Name'] = 'WT'
    df = calc_std(df)

    plates = df['Test plate #'].nunique()

    for plate_no in range(1, plates + 1):
        fig, ax = plt.subplots(figsize=(15,8))
        plate_group = df[df['Test plate #'] == plate_no]
        plate_group = plate_group.dropna(subset=['Avg Exp']).reset_index(drop=True)

        # matplotlib implementation
        # https://stackoverflow.com/questions/58009069/how-to-avoid-overlapping-error-bars-in-matplotlib

        x = plate_group['SampleName'].unique()
        x = [val for val in x for _ in (0, 1)]
        y1 = plate_group.loc[plate_group['Experiment Name'] == 'MT', 'Avg Exp'].tolist()
        yerr1 = plate_group.loc[plate_group['Experiment Name'] == 'MT', 'Exp_std'].tolist()
        y2= plate_group.loc[plate_group['Experiment Name'] == 'WT', 'Avg Exp'].tolist()
        yerr2 = plate_group.loc[plate_group['Experiment Name'] == 'WT', 'Exp_std'].tolist()
        y3= plate_group.loc[plate_group['Experiment Name'] == 'Total', 'Avg Exp'].tolist()
        yerr3 = plate_group.loc[plate_group['Experiment Name'] == 'Total', 'Exp_std'].tolist()

        trans1 = Affine2D().translate(-0.1, 0.0) + ax.transData
        trans3 = Affine2D().translate(+0.1, 0.0) + ax.transData
        er1 = ax.errorbar(x, y1, yerr1, linestyle="none", ecolor='red', transform=trans1)
        er2 = ax.errorbar(x, y2, yerr2, linestyle="none", ecolor='blue')
        er3 = ax.errorbar(x, y3, yerr3, linestyle="none", ecolor='green', transform=trans3)
        ax.scatter(x, y1, c='red', transform=trans1)
        ax.scatter(x, y2, c='blue')
        ax.scatter(x, y3, c='green', transform=trans3)

        plt.xticks(rotation=40, ha='right')
        plt.ylabel('Average Expression')
        plt.title('Plate ' + str(plate_no) + ' Avg Exp')
        plt.suptitle(None)

        plt.grid(True)

        plt.show()


def lm_plot(df):
    # plt.subplots(1, 1, figsize=(15, 10))
    # separate by plates
    df.loc[df['Experiment Name'].str.contains('MT'), 'Experiment Name'] = 'MT'
    df.loc[df['Experiment Name'].str.contains('Total'), 'Experiment Name'] = 'Total'
    df.loc[df['Experiment Name'].str.contains('WT'), 'Experiment Name'] = 'WT'
    df = calc_std(df)

    plates = df['Test plate #'].nunique()
    for plate_no in range(1, plates + 1):
        # plt.figure(figsize=(15,10))
        plate_group = df[df['Test plate #'] == plate_no]
        plate_group = plate_group.dropna(subset=['Avg Exp']).reset_index(drop=True)

        # seaborn implementation
        sns.lmplot("SampleName", "Avg Exp", hue="Experiment Name",
                   data=plate_group, fit_reg=False, x_estimator=np.mean, height=8, aspect=12/8)
        plt.xticks(ha='right')
        plt.ylabel('Average Expression')
        plt.title('Plate ' + str(plate_no) + ' Avg Exp')
        plt.suptitle(None)

        plt.show()


def bar_plot(df):
    # plt.subplots(1, 1, figsize=(15, 10))
    # separate by plates
    df.loc[df['Experiment Name'].str.contains('MT'), 'Experiment Name'] = 'MT'
    df.loc[df['Experiment Name'].str.contains('Total'), 'Experiment Name'] = 'Total'
    df.loc[df['Experiment Name'].str.contains('WT'), 'Experiment Name'] = 'WT'
    df = calc_std(df)

    plates = df['Test plate #'].nunique()
    for plate_no in range(1, plates + 1):
        # fig, ax = plt.subplots()
        # plt.figure(figsize=(15,10))
        plate_group = df[df['Test plate #'] == plate_no]
        plate_group = plate_group.dropna(subset=['Avg Exp']).reset_index(drop=True)

        # sns.barplot(x='SampleName', y='Avg Exp', hue='Experiment Name', data=plate_group)
        g = sns.FacetGrid(data=plate_group, height=8, aspect=12/8)
        g.map(plt.errorbar, 'SampleName', 'Avg Exp', 'Exp_std', fmt='o', elinewidth=1, capsize=5, capthick=1)
        # g.map_dataframe(sns.pointplot, x='SampleName', y='Avg Exp', hue='Experiment Name',palette=sns.color_palette()).add_legend()
        plt.xticks(ha='right')
        plt.ylabel('Average Expression')
        plt.title('Plate ' + str(plate_no) + ' Avg Exp')
        plt.suptitle(None)
        plt.show()
