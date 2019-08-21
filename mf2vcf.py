#! /usr/bin/env python3

import pandas as pd

def read_mf_predictions(mfpred_path):
    mfpred = pd.read_csv(mfpred_path, sep='\t')
    return(mfpred)


def chrom_alt_fields(mfpred):
    '''
    Create #CHROM, POS, REF, and ALT fields from "id"

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing the fields mentioned and, in addition, a
    sample field

    Details:
    There must be one and only one sample name in the sample field otherwise
    ValueError is raised.
    '''
    mat = [y.split(sep='~') for y in mfpred['id']]
    ncol = len(mat[0])
    fields = ['sample', '#CHROM', 'POS', 'REF', 'ALT']
    dtypes = ['category', 'category', 'int64', 'category', 'category']
    if ncol != len(fields):
        raise ValueError('"id" has unexpected number of subfields.')
    # transpose mat to make dictionary from columns
    dic = {fields[i]: [row[i] for row in mat] for i in range(ncol)}
    df = pd.DataFrame(dic)
    d = {y[0]: y[1] for y in zip(fields, dtypes)}
    df = df.astype(d)
    if len(df['sample'].cat.categories) != 1:
        raise ValueError('There must be one and only one sample.')
    def id_field():
        df1 = df[['sample', '#CHROM', 'POS']]
        df2 = pd.DataFrame({'ID': ['MISSING'] * len(df)})
        df3 = df[['REF', 'ALT']]
        res = pd.concat([df1, df2, df3], axis=1)
        return(res)
    df = id_field()
    return(df)


def qual_filter_fields(mfpred):
    '''
    Create QUAL, and FILTER fields from mfpred

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing the fields mentioned
    '''
    F = [y.replace('mosaic', 'PASS') for y in mfpred['prediction']]
    df = pd.DataFrame({'QUAL': 'MISSING', 'FILTER': F})
    return(df.astype('category'))
