import os
import sys

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from modules.utils import merge_redundant_series
from modules.controls import make_dx_disease_gdict

cmap_dir = os.path.dirname(os.path.abspath(__file__))


def make_cmap_genesets(exp, control_exp, valid_genes_file, method='rankdif'):
    with open(valid_genes_file, 'r') as f:
        valid_genes = [line.strip() for line in f]
    up_gsets = make_dx_disease_gdict(exp, control_exp, n_genes=150, up=True, method=method, valid_genes=valid_genes, sets=False)
    dn_gsets = make_dx_disease_gdict(exp, control_exp, n_genes=150, up=False, method=method, valid_genes=valid_genes, sets=False)
    return up_gsets, dn_gsets


def write_cmap_genesets(gs_dicts, out_dir):
    directions = 'up dn'.split()
    for i, gs_dict in enumerate(gs_dicts):
        for sample, genes in gs_dict.items():
            with open(os.path.join(out_dir, '{}.{}.txt'.format(sample, directions[i])), 'w') as f:
                f.write('\n'.join(genes))


def read_cmap_gct(gct_path, score_only=True, pert_type='CP'):
    results = pd.read_table(gct_path, skiprows=2, index_col=0).T
    new_index = list(results.cell_id)
    new_index[:7] = 'pert_type pert_name subclasses targets classes pc ts_pc'.split()
    results.index = new_index
    results.drop('cell_id _id name pc pc_selection'.split(), axis=1, inplace=True)
    if pert_type is not None:
        results = results.T
        results = results[results.pert_type == pert_type].T
    if score_only:
        results.columns = results.loc['pert_name']
        results.drop(['pert_type',
                      'pert_name',
                      'subclasses',
                      'targets',
                      'classes',
                      'pc',
                      'ts_pc',
                      'PC3',
                      'VCAP',
                      'A375',
                      'A549',
                      'HA1E',
                      'HCC515',
                      'HT29',
                      'MCF7',
                      'HEPG2'], inplace=True)
        results.columns.name = 'drug_name'
        results.index = ['score']
        results = results.T
        results.loc[:, 'score'] = results.loc[:, 'score'].astype(float)
        results = merge_redundant_series(pd.DataFrame(results), axis=0)
        results.sort_values(by='score', inplace=True)

    return results * -1 / 100  # so that higher is better. Otherwise, anticorrelated drugs are better
