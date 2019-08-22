#! /usr/bin/env python3

__format_spec__ = '.4g'

__info__ = {\
        'mappability':
        {'ID': 'MAPPABILITY', 'Number': 1, 'Type': 'Float', 'Description': 'NA'},
        'type':
        {'ID': 'TYPE', 'Number': 1, 'Type': 'String', 'Description': 'NA'},
        'length':
        {'ID': 'LENGTH', 'Number': 1, 'Type': 'Integer', 'Description': 'NA'},
        'GCcontent':
        {'ID': 'GCCONTENT', 'Number': 1, 'Type': 'Float', 'Description': 'NA'},
        'context':
        {'ID': 'CONTEXT', 'Number': 1, 'Type': 'String', 'Description': 'NA'}}

__format__ = {\
        'AF':
        {'ID': 'AF', 'Number': 1, 'Type': 'Float', 'Description': 'NA'},
        'dp':
        {'ID': 'DP', 'Number': 1, 'Type': 'Integer', 'Description': 'NA'}}

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
    FILTER = [y.replace('mosaic', 'PASS') for y in mfpred['prediction']]
    df = pd.DataFrame({'QUAL': 'MISSING', 'FILTER': FILTER})
    return(df.astype('category'))


def info_fields(mfpred):
    '''
    Create the INFO fields from mfpred

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing the INFO fields separated by a ";"

    Details:
    The __info__  variable (a dict) is used to define the ID, Number, Type and
    Description of a particular INFO field.  The __format_spec__ variable is
    used to format all INFO fields of Type "float".
    '''
    def one_info_field(mfid='mappability'):
        ID = __info__[mfid]['ID']
        Type = __info__[mfid]['Type']
        def info_float(x):
            val = ID + '=' + str(format(x, __format_spec__))
            return(val)
        def info_integer(x):
            val = ID + '=' + str(x)
            return(val)
        def info_string(x):
            val = ID + '=' + x
            return(val)
        if Type == 'Float':
            helper = info_float
        elif Type == 'Integer':
            helper = info_integer
        elif Type == 'String':
            helper = info_string
        elif Type == 'Flag':
            raise ValueError('Flag INFO field type is not yet implemented')
        data = [helper(y) for y in mfpred[mfid]]
        return(data)
    fields = __info__.keys()
    val = [one_info_field(mfid=f) for f in fields]
    val = [[row[i] for row in val] for i in range(len(fields))]
    val = [';'.join(y) for y in val]
    df = pd.DataFrame({'INFO': val})
    return(df)


def format_sample_fields(mfpred):
    def one_sample_field(mfpred, mfid='AF'):
        ID = __format__[mfid]['ID']
        Type = __format__[mfid]['Type']
        def sample_float(x):
            return(str(format(x, __format_spec__)))
        def sample_general(x):
            return(str(x))
        if Type == 'Float':
            helper = sample_float
        else:
            helper = sample_general
        val = [helper(y) for y in mfpred[mfid]]
        return(val)
    fields = __format__.keys()
    ids = [__format__[f]['ID'] for f in fields]
    FORMAT = ':'.join(ids)
    vals = [one_sample_field(mfpred, mfid=f) for f in fields]
    vals = [':'.join(y) for y in zip(*vals)]
    df = pd.DataFrame({'FORMAT': FORMAT, 'SAMPLE': vals})
    return(df)



def data_lines(mfpred):
    '''
    Creates all data lines (i.e. the body) of the VCF

    Parameter:
    mfpred: the pandas DataFrame from read_mf_predictions

    Returns:
    a pandas DataFrame containing all fields
    '''
    funs = [chrom_alt_fields, qual_filter_fields, info_fields, format_sample_fields]
    l = [f(mfpred) for f in funs]
    df = pd.concat(l, axis=1)
    sample = df['sample'][0]
    df1 = df[df.columns[1:-1]]
    df2 = pd.DataFrame({sample: df['SAMPLE']}, index=range(len(df)))
    val = pd.concat([df1, df2], axis=1)
    return(val)


