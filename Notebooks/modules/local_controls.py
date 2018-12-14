import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from modules.utils import reduce_to_common_idxs, rank_normalize

def make_dx_disease_gdict(exp_df, control_exp_df, n_genes=150, up=True, method='rankdif', valid_genes=None, sets=True):
    disease_gdict = {}
    control_exp_olap, exp_olap = reduce_to_common_idxs([control_exp_df, exp_df])
    if method == 'rankdif':
        control_ranks, ranks = [rank_normalize(df) for df in [control_exp_olap, exp_olap]]
        rankby_df = ranks.subtract(control_ranks.iloc[0])
    elif method == 'logfold':
        ratios = exp_df.divide(control_exp_df.iloc[0])
        rankby_df = pd.DataFrame(np.log2(ratios.values), index=ratios.index, columns=ratios.columns)
    else:
        raise ValueError('method "{}" not supported; try one of ["rankdif", "logfold"]'.format(method))
    if valid_genes is None:
        valid_genes_in_exp = exp_olap.columns.tolist()
    else:
        valid_genes_in_exp = []
        for g in exp_olap.columns:
            for subg in g.split(' /// '):
                if subg in valid_genes:
                    valid_genes_in_exp.append(subg)
                    break
    for sample in exp_df.index:
        topn = rankby_df.loc[sample, valid_genes_in_exp].sort_values(ascending=False if up else True).head(n_genes).index.tolist()
        if sets:
            topn = set(topn)
        disease_gdict[sample] = topn

    return disease_gdict


def load_control_exp(exp_drug_suggestion_controls_dir, control='neural_stem'):
    control_exp_file = os.path.join(exp_drug_suggestion_controls_dir, '{}.csv'.format(control))
    control_exp = pd.read_csv(control_exp_file, index_col=0)
    return control_exp
