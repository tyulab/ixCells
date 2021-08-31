import click
import os
import pandas as pd


def file_split(file):
    s = file.split('.')
    name = '.'.join(s[:-1])  # get directory name
    return name


def getsheets(inputfile, fileformat):
    name = file_split(inputfile)
    try:
        os.makedirs(name)
    except:
        pass

    df1 = pd.ExcelFile(inputfile)
    for x in df1.sheet_names:
        print(x + '.' + fileformat, 'Done!')
        df2 = pd.read_excel(inputfile, sheet_name=x)
        filename = os.path.join(name, x + '.' + fileformat)
        if fileformat == 'csv':
            df2.to_csv(filename, index=False)
        else:
            df2.to_excel(filename, index=False)
    print('\nAll Done!')


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('inputfile')
@click.option('-f', '--format', type=click.Choice([
    'xlsx', 'csv']), default='xlsx', help='Default xlsx.')
def cli(inputfile, format):
    '''Convert a Excel file with multiple sheets to several file with one sheet.
    Examples:
    \b
        getsheets filename
    \b
        getsheets filename -f csv
    '''
    if format == 'csv':
        getsheets(inputfile, 'csv')
    else:
        getsheets(inputfile, 'xlsx')


cli()