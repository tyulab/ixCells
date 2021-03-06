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
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
import matplotlib.backends.backend_pdf

from utils.stats import *
from utils.functions import *
import config

def type_hist(df, col):
    pdf = matplotlib.backends.backend_pdf.PdfPages(get_folder() + 'Z Score Range Plots/Exp zscore hist.pdf')

    # hide naive/control
    df = hide_naive_control(df)

    # separate into WT, MT, total
    df_no_na = df.dropna(subset=[col]).reset_index(drop=True)
    types = ['WT','MT','Total']


    for i in range(len(types)):
        type_df = df_no_na[df_no_na['Experiment Name'].str.contains(types[i], case=False)]
        ax = type_df.hist(column=col, bins=20).ravel()[0].get_figure()
        plt.title("Exp zscore ranges for "+types[i])
        plt.xlabel('Exp zscore range')
        # save the dfs in case
        type_df.to_csv(get_folder()+"Z Score Range Plots/Exp zscore "+types[i]+".csv", index=False)
        if config.SAVE_FIG_AS == 'pdf':
            # fig = ax[0].get_figure()
            pdf.savefig(ax)
        else:
            plt.savefig(get_folder() +'Z Score Range Plots/' + 'Exp zscore hist ' + types[i] + '.png')
        if config.SHOW_PLOTS:
            plt.show()

    pdf.close()

# create histogram for controls with matplotlib
def control_hist(df):
    # pos = "Ionis1375651"
    # neg = "Ionis676630"
    # separate _3 from _10
    df = df[['Experiment Name', 'Position', 'SampleName', 'ASO Microsynth ID', 'Exp', 'Test plate #']]# .dropna(subset=['Exp'])
    controls = df['SampleName'].str.contains('Ionis|Naive')
    df[controls].to_csv(get_folder()+"Control Plots/controls.csv", index=False)
    control_10 = df[controls & df['SampleName'].str.contains('_10')].reset_index(drop=True)
    control_3 = df[controls & df['SampleName'].str.contains('_3')].reset_index(drop=True)
    # print(control_10)
    control_10.to_csv(get_folder()+"Control Plots/controls_10.csv", index=False)
    control_3.to_csv(get_folder()+"Control Plots/controls_3.csv", index=False)
    if not control_10.empty:
        hist_plates(control_10, "hist_10", title='10mmol Plates')
    if not control_3.empty:
        hist_plates(control_3, "hist_3", title='3mmol Plates')

# helper function to create histograms by plates
def hist_plates(control, file_prefix="hist", title='Plates'):
    pdf = matplotlib.backends.backend_pdf.PdfPages(get_folder() + 'Control Plots/'+file_prefix + '.pdf')

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

    plates = control['Test plate #'].unique()
    plates.sort()
    plate_count = len(plates)
    step = 9 # min(plate_count, 9)
    layout = (3,3)

    for i in tqdm(range(0,plate_count,step)):
        # fig, ax1 = plt.subplots(nrows=layout[0], ncols=layout[1])
        plt.margins(x=0, y=0)
        rem = plate_count-i
        # TODO: change between somehow?
        # axes = pos_df[pos_df['Test plate #'].between(i, i+step-1)].hist( \
        axes = pos_df[pos_df['Test plate #'].isin(plates[i:i+step-1])].hist( \
            column='Exp',by=pos_df['Test plate #'], layout=layout, # ax = ax1,\
             sharex=True, sharey=True, figsize=(6,7), xrot=0, \
            alpha=0.5, label='Ionis1375651', bins=5)
        axes = axes.ravel()[:min(rem,step)]

        # plt.yticks(np.linspace(0, 10, num=6, endpoint=True))
        plt.ylim(bottom=0,top=10)

        if neg_df is not None:
            plt.xlim(left=0.0, right=1.75)
            plt.xticks(np.linspace(0,2,num=5,endpoint=True))
            # neg_df[neg_df['Test plate #'].between(i, i+step-1)].hist( \
            neg_df[neg_df['Test plate #'].isin(plates[i:i + step - 1])].hist( \
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

        # https://stackoverflow.com/a/67202037
        fig = axes.ravel()[0].get_figure()
        if config.SAVE_FIG_AS == 'pdf':
            pdf.savefig(fig)
        else:
            plt.savefig(get_folder()+'Control Plots/'+file_prefix+' '+str(i)+'-'+str(i+min(rem,step)-1)+'.png')
    if config.SHOW_PLOTS:
        plt.show()
    pdf.close()

# helper function to create scatter plot for the tiers df
def tier_scatter(df, title='Tiers by MT/WT'):
    # Create nums for MT and WT
    df['MT Avg Exp'] = df['MT Avg Exp']*100
    df['WT Avg Exp'] = df['WT Avg Exp']*100

    # cmap = colors.ListedColormap(['r', 'g', 'b', 'c'])
    # bounds = [0, 10, 20]
    # norm = colors.BoundaryNorm(bounds, cmap.N)

    df.plot.scatter(x='WT Avg Exp', y = 'MT Avg Exp') # alpha=0.5

    plt.title(title)
    plt.savefig(get_folder()+'Tiers/scatter_uncropped.pdf')

    df = df[df['MT Avg Exp'] <= 110]
    df = df[df['WT Avg Exp'] <= 110]
    # ax1 = df.plot.scatter(x='WT Avg Exp', y = 'MT Avg Exp', c=df['MT Avg Exp'], cmap=cmap)
    plt.xlim(left=0.0, right=100)
    plt.ylim(bottom=0, top=100)
    plt.xticks(np.linspace(0,110,num=12,endpoint=True))
    plt.yticks(np.linspace(0,110,num=12,endpoint=True))
    plt.grid(True)

    plt.savefig(get_folder()+'Tiers/scatter.pdf')
    if config.SHOW_PLOTS:
        plt.show()

# create r2 plots
def r2_plot(y_test,y_predicted, title, output_file):
    fig, ax = plt.subplots()
    ax.scatter(y_test, y_predicted)
    # ax.plot([y_test.min(), y_test.max()], [y_predicted.min(), y_predicted.max()], 'k--', lw=4)
    ax.set_xlabel('B1')
    ax.set_ylabel('B2')
    # regression line
    y_test, y_predicted = y_test.reshape(-1, 1), y_predicted.reshape(-1, 1)
    model= LinearRegression().fit(y_test, y_predicted) # Added by AS
    ax.plot(y_test, model.predict(y_test), color='red', linewidth=0.75) # Added by AS
    r2 = str(model.score(y_test, y_predicted)) # Added by AS
    ax.set_title(title + ' R2: ' + r2) # Added by AS
    #ax.plot(y_test, LinearRegression().fit(y_test, y_predicted).predict(y_test), color='red', linewidth=0.75) #KD
    #ax.set_title(title+' R2: ' + str(r2_score(y_test, y_predicted))) #KD
    #plt.savefig(get_folder()+title+'_R2.pdf')
    # TODO: new file for each round
    output_file.write("{}\t{}\n".format(title, r2))
    if config.SHOW_PLOTS:
        plt.show()
    plt.close()
    return fig

# create error bar plots similar to prism output
def error_plot(df):
    spacing = 0.1

    # plt.subplots(1, 1, figsize=(15, 10))
    # keep only _10 mmol
    # separate by plates
    df.loc[df['Experiment Name'].str.contains('MT'), 'Experiment Name'] = 'MT'
    df.loc[df['Experiment Name'].str.contains('Total'), 'Experiment Name'] = 'Total'
    df.loc[df['Experiment Name'].str.contains('WT'), 'Experiment Name'] = 'WT'
    df = calc_std(df)
    # # SEM
    # df = calc_sem(df)

    # TODO: iterate through plates as string
    plates = df['Test plate #'].unique()

    pdf_name = get_folder()+config.ROUND+' kcnq2_exp'
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name+'.pdf')
    # if plates > 20:
    #     pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name+'_P1-20.pdf')

    for plate_no in tqdm(range(len(plates))):
        fig, ax = plt.subplots(figsize=(15,8))

        # TODO: make sure no plates are skipped (Round 1 issue)
        plate = df[df['Test plate #'] == plates[plate_no]]
        plate = plate.dropna(subset=['Avg Exp']).reset_index(drop=True)
        # plate_group = plate

        # re-order naive to left of graph?
        naive = plate['ASO Microsynth ID'].str.contains('Naive', case=False)
        naive_index = plate.index[naive].tolist()
        # re-order ionis
        ionis = plate['ASO Microsynth ID'].str.contains('Ionis', case=False)
        ionis_index = plate.index[ionis].tolist()
        other_index = plate.index[~naive & ~ionis].tolist()
        # plate = pd.concat([plate[naive], df[~naive]], ignore_index=True).copy()
        plate_group = plate.reindex(naive_index + ionis_index + other_index).reset_index(drop=True).copy()


        # matplotlib implementation
        # https://stackoverflow.com/questions/58009069/how-to-avoid-overlapping-error-bars-in-matplotlib
        # TODO: delete some unused list
        y = []
        y_mean = []
        yerr = []
        err = []
        scatter = []
        colors=['r','b','g']
        types=['MT', 'WT', 'Total']

        for idx, t in enumerate(types):
            # unique names of samples
            x = plate_group[plate_group['Experiment Name'].str.contains(t, case=False)]['ASO Microsynth ID'].unique()
            # all occurences of avg exp points
            x2 = plate_group[plate_group['Experiment Name'].str.contains(t, case=False)]['ASO Microsynth ID'].tolist()
            y.append(plate_group[plate_group['Experiment Name'].str.contains(t, case=False)]['Avg Exp'].tolist())
            y_mean.append(plate_group[plate_group['Experiment Name'].str.contains(t, case=False)].groupby('ASO Microsynth ID',sort=False)['Avg Exp'].mean().tolist())
            yerr.append(plate_group[plate_group['Experiment Name'].str.contains(t, case=False)].groupby('ASO Microsynth ID',sort=False)['Exp_std'].first().tolist())
            transform = Affine2D().translate(spacing*idx-spacing, 0.0) + ax.transData

            scatter.append(ax.scatter(x2, y[idx], c=colors[idx], s=25, transform=transform))
            err.append(ax.errorbar(x, y_mean[idx], yerr[idx], linestyle='none', c=colors[idx],fmt='s', markersize=2.5, transform=transform ))


        # TODO: change alpha value for dosage
        # bool_10 = df['SampleName'].str.contains('_10').tolist()
        # alphas = [1.0 if x else 0.5 for x in bool_10]

        plt.xticks(rotation=40, ha='right')
        plt.grid(True, axis='y', linestyle=':', linewidth=2, zorder=32)
        plt.ylabel('KCNQ2 Expression')
        if config.RENORMALIZE:
            plt.ylim(0, 3)
            plt.yticks(np.linspace(0, 3, 7))
        else:
            plt.ylim(0,1.5)
            plt.yticks(np.linspace(0,1.5,7))
        # plt.xlabel('ASOs')
        # if renormalized to neg change title
        if config.RENORMALIZE:
            title = config.ROUND + ' Plate ' + str(plates[plate_no]) + ' Normalized to Negative Control'
        else:
            title = config.ROUND + ' Plate ' + str(plates[plate_no]) + ' Normalized to Naive'
        plt.title(title)
        plt.suptitle(None)
        plt.tight_layout()
        plt.legend(scatter,types)
        if config.SHOW_PLOTS:
            plt.show()

        pdf.savefig(fig)

        plt.close()

        # if plate_no % 20 == 0:
        #     pdf.close()
        #     pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name+'_P'+str(plate_no+1)+'-'+str(min(plate_no+19,plates))+'.pdf')

    pdf.close()