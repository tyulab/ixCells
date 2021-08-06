#!/usr/bin/env python
from utils.functions import *
from utils.plots import *
from utils.stats import *
import config

# limit on z score range to filter out before making list of tiers
# threshold for output
# normalize to

# make output file
# steps- renormalize ddCt, recalculate Exp, calculate abs z scores, drop outliers, avg z scores on samples, assign plates, avg on plates
def create_output():
    # specify path to folder containing all csvs and plate sheet + output to be ignored
    # TODO: specify round in config file
    path = config.ROUND+"\data"
    # For Round 1
    plate_file = "ixCells_Round 1_2021-06-22_TN09_551ASOs_plate id adjusted.csv"
    output_file = get_folder() + "/output.csv"
    csv = get_csv(path, plate_file)

    pd.set_option('display.max_columns', None)

    df_list = []
    print('found ' + str(len(csv)) + ' csv files')
    # read all csvs into dataframes and concatenate
    for f in tqdm(csv):
        if plate_file not in f:
            df = pd.read_csv(f, encoding='latin-1')
            df.columns = df.columns.str.strip()
            # rename short description to microsynth id to have consistent columns, original naive normalized cols
            df = df.rename(columns={'ASO Short Description': 'ASO Microsynth ID'})
            # drop fully empty rows
            df = df.dropna(how='all').reset_index(drop=True)
            # redo dCt from Crossing Point
            df = redo_dct(df)
            # drop irrelevant columns
            df = df[['Experiment Name', 'Position', 'SampleName', 'ASO Microsynth ID', 'dCt', 'ddCt', 'Exp']]

            # if want to renormalize
            if config.RENORMALIZE:
                # rename old columns and get negative control ddct
                df = df.rename(columns={'ddCt': 'ddCt from Naive', 'Exp':'Exp from Naive'})
                neg_ddct(df)

            # rename naive
            df = df.replace(to_replace='^Na.*e', value='Naive', regex=True)


            df_list.append(df)

    df = pd.concat(df_list, ignore_index=True)
    # pd.set_option('display.max_rows', df.shape[0]+1)

    # Calculate exp
    calc_exp(df)
    # renormalizing
    if config.RENORMALIZE:
        calc_exp(df, col='ddCt from Naive', output_col='Exp from Naive')

    # Calculate z scores for Exp
    abs_zscore(df)
    # drop outliers
    # set threshold on z exp zscores
    threshold = flag_outliers(df, col='Exp_zscore', greater_than=config.Z_SCORE_THRESHOLD)
    dropped = df[threshold].reset_index(drop=True)
    dropped.to_csv(get_folder()+"/output_dropped.csv", index=False)
    df[~threshold].reset_index(drop=True).to_csv(get_folder()+"/output_filtered.csv", index=False)
    df.loc[threshold,'Exp'] = np.nan
    # df = df[~threshold].reset_index(drop=True)
    # Average z scores
    # todo: remove
    # avg_zscore(df)
    assign_plates(df, plate_file)
    avg_plates(df) # can separate this part later?
    # print(df)
    # export to output file
    df.to_csv(output_file, index=False)

# create avg exp file
def create_avg_exp():
    # default path for table to read from
    table_file = get_folder()+"/output.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    df = df[['Experiment Name','SampleName','ASO Microsynth ID','Exp_zscore','Exp','Test plate #']]
    pd.set_option('display.max_columns', None)

    # # hide naive / controls from avg exp
    # df = hide_naive_control(df)

    # calculate avg exp, range
    df = exp_zscore_range(df)
    df.to_csv(get_folder()+"/avg_exp.csv", index=False)
    # make histogram from ranges
    type_hist(df, 'Avg Exp_zscore range')

def create_tiers():
    # default path for table to read from
    table_file = get_folder()+"/avg_exp.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # hide naive / controls from avg exp
    df = hide_naive_control(df)

    # remove 3mmol and total from view
    df = remove_3(df)
    df = remove_total(df)

    # flag ranges > threshold and drop
    df = df.dropna(subset=['Avg Exp_zscore range']).reset_index(drop=True)
    ranges = flag_outliers(df, col='Avg Exp_zscore range', greater_than=config.TIERS_THRESHOLD)
    dropped = df[ranges].sort_values('Avg Exp_zscore range', ascending=True).reset_index(drop=True)
    dropped_groups = dropped.groupby('SampleName').first().sort_values(['SampleName'], na_position='last')
    dropped_groups = dropped_groups[['ASO Microsynth ID', 'Test plate #']]
    dropped_groups['Threshold: ' + str(config.TIERS_THRESHOLD)] = np.nan
    dropped_groups.to_csv(get_folder()+"Tiers/tiers_dropped.csv", index=True)

    # remove anything in dropped samples from tiering
    dropped_samples = dropped['SampleName'].unique()
    df = df[~df['SampleName'].isin(dropped_samples)].reset_index(drop=True)

    tiers = tierlist(df)

    tiers.sort_values(['Tier','SampleName'], ascending=True, inplace=True)
    tiers.to_csv(get_folder()+"Tiers/tiers.csv", index=True)


# make control histograms
def create_control_hist():
    # default path for table to read from
    table_file = get_folder()+"output.csv"
    plate_output = "plates_table.csv"
    # store in df
    df = pd.read_csv(table_file, encoding='latin-1')
    pd.set_option('display.max_columns', None)

    # TODO: flag any controls across plates > 1.55
    print('control plots:')
    outliers = flag_outliers(df, col='Exp', greater_than=1.55)
    df = df[~outliers].reset_index(drop=True)
    # histogram
    control_hist(df)


def create_tier_plots():
    print("create tier plots")
    file = get_folder()+"Tiers/tiers.csv"
    # lst of column names which needs to be string
    lst_str_cols = ['Tier']
    # use dictionary comprehension to make dict of dtypes
    dict_dtypes = {x: 'str' for x in lst_str_cols}
    # use dict on dtypes
    df = pd.read_csv(file, dtype=dict_dtypes)
    # scatter plot
    tier_scatter(df)

# r squared analysis
def create_r_squared():
    file = get_folder()+"avg_exp.csv"
    # store in df
    df = pd.read_csv(file, encoding='latin-1')
    pd.set_option('display.max_columns', None)
    print('create r squared plots')

    # hide naive / controls from avg exp
    df = hide_naive_control(df)

    # keep only 10mmol
    df = df[df['SampleName'].str.contains('_10')]

    pdf = matplotlib.backends.backend_pdf.PdfPages(get_folder() + 'r_squared.pdf')
    # get replicates
    for type in ['WT','MT']:
        type_df = get_type(df, type=type)
        b1,b2 = get_replicates(type_df)
        fig = r2_plot(b1,b2,title=type)
        pdf.savefig(fig)

    pdf.close()

def create_error_bars():
    file = get_folder()+"avg_exp.csv"
    # store in df
    df = pd.read_csv(file, encoding='latin-1')
    pd.set_option('display.max_columns', None)
    print('create error bar plots')
    error_plot(df)