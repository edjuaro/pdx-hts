from functools import reduce
from collections import defaultdict, Counter
import itertools

import pandas as pd
import numpy as np
from pandas import ExcelWriter
import networkx as nx
from scipy.stats import rankdata
from sklearn.preprocessing import scale


def make_col2col_dict(tabular_file, from_col, to_col, sep='\t'):
    with open(tabular_file, 'r') as f:
        cols = f.readline().strip().split(sep)
    for col in [from_col, to_col]:
        if col not in cols:
            raise Exception("{} not in columns".format(col))
    df = pd.read_table(tabular_file, sep=sep, usecols=[from_col, to_col], dtype=str).dropna(axis=0)
    col2col = dict(zip(df[from_col], df[to_col]))
    return col2col


def multi_intersection(iterables):
    return reduce(set.intersection, map(set, iterables))


def multi_union(iterables):
    return reduce(set.union, map(set, iterables))


def save_xls(list_dfs, sheetnames, xls_path, fit_col_width=True, max_col_width=30):
    writer = ExcelWriter(xls_path)
    for n, df in enumerate(list_dfs):
        df.to_excel(writer, sheetnames[n])
        if fit_col_width:
            worksheet = writer.sheets[sheetnames[n]]  # pull worksheet object
            index = df.index
            index_name = index.name if index.name is not None else ''
            max_len = max([len(x) for x in index] + [len(index_name)]) + 1
            worksheet.set_column(0, 0, min(max_len, max_col_width))
            for idx, col in enumerate(df):  # loop through all columns
                series = df[col]
                max_len = max((
                    series.astype(str).map(len).max(),  # len of largest item
                    len(str(series.name))  # len of column name/header
                )) + 1  # adding a little extra space
                worksheet.set_column(idx + 1, idx + 1, min(max_len, max_col_width))
    writer.save()


def reverse_item2set_dict(p2g):
    g2p = defaultdict(set)
    for p, gs in p2g.items():
        for g in gs:
            g2p[g].add(p)
    return g2p


def standardize_idxs_and_dedup(exp, dr):
    # exp.index = [idx.split('_')[0] for idx in exp.index]
    exp, dr = reduce_to_common_idxs([exp, dr], axis=0)
    dr = merge_redundant_series(dr, axis=0)
    dr = merge_redundant_series(dr, axis=1)
    exp = merge_redundant_series(exp, axis=0)
    exp = merge_redundant_series(exp, axis=1)
    return exp, dr


def expand_cols_to_synonymns(dfs, syn_delim=' /// '):
    dfcs = [df.copy() for df in dfs]
    syn_graph = nx.Graph()
    for df in dfcs:
        for col in df.columns:
            syns = col.split(syn_delim)
            syn_graph.add_nodes_from(syns)
            for s1, s2 in itertools.combinations(syns, 2):
                syn_graph.add_edge(s1, s2)
    syn2syns = {}
    ccs = nx.connected_components(syn_graph)
    for cc in ccs:
        syn_string = syn_delim.join(sorted(cc))
        for syn in cc:
            syn2syns[syn] = syn_string
    for i, df in enumerate(dfcs):
        newcols = []
        for col in df.columns:
            syn1 = col.split(syn_delim)[0]
            newcols.append(syn2syns[syn1])
        df.columns = sorted(newcols)
        dfcs[i] = merge_redundant_series(df, axis=1)
    return dfcs


def reduce_to_common_idxs(dfs, axis=1, syns=True):
    dfcs = [df.T.copy() if axis == 0 else df.copy() for df in dfs]
    if syns:
        dfcs = expand_cols_to_synonymns(dfcs)
    idxs = [set(df.columns) for df in dfcs]
    idxs_in_common = sorted(set.intersection(*idxs))
    dfcs = [df.loc[:, idxs_in_common] for df in dfcs]
    return [df.T if axis == 0 else df for df in dfcs]


def merge_redundant_series(orig_df, axis=0, method='mean', mask_nans=True):
    df = orig_df.T if axis == 1 else orig_df
    idx_counts = Counter(df.index)
    already_single = []
    mult_idxs = []
    for idx, count in idx_counts.items():
        (already_single if count == 1 else mult_idxs).append(idx)
    single_subdf = df.loc[already_single]
    if len(mult_idxs) == 0:
        new_df = single_subdf
    else:

        series_list = []
        for uid in mult_idxs:
            subdf = df.loc[uid]
            if isinstance(subdf, pd.DataFrame):
                if mask_nans:
                    mask = np.isnan(subdf)
                    vals = np.ma.masked_array(subdf, mask)
                else:
                    vals = subdf.values
                if method == 'mean':
                    series = vals.mean(axis=0)
                elif method == 'max':
                    series = vals.max(axis=0)
                elif method == 'min':
                    series = vals.min(axis=0)
                elif method == 'first':
                    series = vals[0]
                elif method == 'sum':
                    series = vals.sum(axis=0)
                else:
                    raise ValueError("{} not a supported method for merge".format(method))
                series = pd.Series(series, index=subdf.columns)
                series.name = uid
            else:
                series = subdf
            series_list.append(series)

        dedup_df = pd.concat(series_list, axis=1).T
        # print(single_subdf.shape, dedup_df.shape)
        new_df = pd.concat([dedup_df, single_subdf], axis=0)
    return new_df.T if axis == 1 else new_df


def read_gct(gct_path, use_description=False, description_func=lambda x: x, transpose=True, dropna_axes=[]):
    df = pd.read_table(gct_path, skiprows=2, index_col=0)
    for axis in dropna_axes:
        df = df.dropna(axis=axis)
    if 'description' in df.columns:
        dcol = 'description'
    elif 'Description' in df.columns:
        dcol = 'Description'
    if use_description:
        df.index = [description_func(d) for d in df[dcol]]
    df = df.drop(dcol, axis=1)
    df = df.dropna(axis=0)
    if transpose:
        df = df.T
    return df


def write_gct(tdf, path):
    df = tdf.T
    n_rows, n_cols = df.shape
    with open(path, 'w') as f:
        f.write("#1.2\n")
        f.write("{}\t{}\n".format(n_rows, n_cols))
        header_row = 'Name\tDescription\t' + '\t'.join(df.columns) + '\n'
        f.write(header_row)
        for i in range(n_rows):
            name = df.index[i]
            values = df.iloc[i]
            row = '{}\t{}\t{}\n'.format(name, name, '\t'.join(map(str, values)))
            f.write(row)


def all_unique_values(set_dict):
    return set(itertools.chain.from_iterable(set_dict.values()))


def rank_normalize(x_orig, ascending=True, method='dense', norm_to_max=True, add_jitter=False,
                   jitter_scale=1E-9):  # , ties=False
    # Todo: profile two methods, see if rankdata is slower, if so efficiently check for ties within each row
    x = x_orig.copy()
    if add_jitter:
        jitter = jitter_scale * np.random.uniform(size=x.shape)
        x += jitter
    n, m = x.shape
    s = x.values if isinstance(x, pd.DataFrame) else x
    ranks = np.array([rankdata(s[i], method=method) for i in range(n)]) - 1
    if not ascending:
        ranks = (m - ranks) - 1
    if norm_to_max:
        ranks = ranks / (m - 1)
    return pd.DataFrame(ranks, index=x.index, columns=x.columns) if isinstance(x, pd.DataFrame) else ranks


def quantile_normalize(x):
    """
    # Wikipedia example:
    x = np.array([[5,    4,    3],
    [2,    1,    4],
    [3,    4,    6],
    [4,    2,    8]])
    index = list('ABCD')
    columns = range(1,4)
    x = pd.DataFrame(x, index=index, columns=columns).T
    """
    n, m = x.shape
    s = x.values if isinstance(x, pd.DataFrame) else x
    ranks = np.array([rankdata(s[i], method='dense') for i in range(n)]) - 1
    sorted_s = np.sort(s, axis=1)
    col_means = np.mean(sorted_s, axis=0)
    ranks_copy = ranks.copy().astype(float)
    for rank in range(ranks.shape[1]):
        ranks_copy[ranks == rank] = col_means[rank]
    return pd.DataFrame(ranks_copy, index=x.index, columns=x.columns) if isinstance(x, pd.DataFrame) else ranks_copy


def scale_df(df):
    scaled_values = scale(df.values)
    return pd.DataFrame(scaled_values, index=df.index, columns=df.columns)


def scale_features_between_zero_and_one(df, fillna=None):
    dfc = df.copy()
    dfc = dfc.subtract(df.min(axis=0))
    mx = dfc.max(axis=0)
    dfc = dfc.divide(mx)
    if fillna is not None:
        dfc = dfc.fillna(fillna)
    return dfc


def permute_columns(x, rs=None):
    if rs is None:
        rs = np.random.RandomState()
    ix_i = rs.random_sample(x.shape).argsort(axis=0)
    ix_j = np.tile(np.arange(x.shape[1]), (x.shape[0], 1))
    return x[ix_i, ix_j]
