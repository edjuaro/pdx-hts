import os
import sys
from collections import defaultdict, Counter

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from modules.utils import merge_redundant_series, expand_cols_to_synonymns, all_unique_values


def to_cdfs_df(df, pad=1):
    vals = df.values.ravel()
    cdfs = to_cdfs(vals, pad=pad).reshape(df.shape)
    cdf_df = pd.DataFrame(cdfs, index=df.index, columns=df.columns)
    return cdf_df


def to_cdfs(x, pad=1):
    a = np.array(x)
    nonnan = a[~np.isnan(a)]
    n_nonnan = len(nonnan)
    val2idxs = defaultdict(list)
    for i, val in enumerate(x):
        val2idxs[val].append(i)
    counts = Counter(nonnan)
    cdfs = np.empty(len(x))
    cdfs.fill(np.nan)
    n_le = 0
    for val in sorted(counts.keys()):
        count = counts[val]
        n_le += count
        idxs = val2idxs[val]
        cdfs[idxs] = n_le
    return cdfs / (n_nonnan + pad)


def log_oddsify(df, pad=1):
    cdfs_df = to_cdfs_df(df, pad=pad)
    probs = cdfs_df.values
    log_odds_df = pd.DataFrame(np.log(np.divide(probs, 1 - probs)), index=df.index, columns=df.columns)
    return log_odds_df


def make_union_drug_matrices(dfs, n2cidss, merge_methods=None):
    if merge_methods is None:
        merge_methods = ['mean'] * len(dfs)
    dfcs = [df.copy() for df in dfs]
    for df, n2cids in zip(dfcs, n2cidss):
        newcols = []
        for col in df.columns:
            cids = n2cids[col] if col in n2cids else []
            newcols.append(' /// '.join(map(str, [col] + cids)))
        df.columns = newcols
    dfcrs = expand_cols_to_synonymns(dfcs)
    nodups = [merge_redundant_series(df, axis=1, method=merge_methods[i]) for i, df in enumerate(dfcrs)]
    return nodups


def simplify_syn_index(syn_index):
    new_idx = []
    for idx in syn_index:
        syns = idx.split(' /// ')
        min_len = 10000
        best_syn = None
        for syn in syns:
            if len(syn) > 0 and not syn.isdigit() and not syn.startswith('SBI'):
                synlen = len(syn)
                if synlen < min_len:
                    min_len = synlen
                    best_syn = syn
        if best_syn is None:
            for syn in syns:
                if syn.isdigit():
                    best_syn = syn
        if best_syn is None:
            best_syn = np.array(syns)[np.argmin([len(syn) for syn in syns])]
        if '_' in best_syn:
            best_syn = best_syn.split('_')[1]
        new_idx.append(best_syn)
    return new_idx


def load_reasonable_drugs(drug_annotation_dir):
    reasonable_drugs_file = os.path.join(drug_annotation_dir, 'clinically_relevant_drugs.csv')
    with open(reasonable_drugs_file, 'r') as f:
        rdrugs = [row.strip() for row in f.readlines()]
    return rdrugs




def select_reasonable_drugs(other_drug_to_cids, rdrug_to_cids):
    reasonable_cids = all_unique_values(rdrug_to_cids)
    cid_to_other_drugs = defaultdict(set)
    for other_drug, cids in other_drug_to_cids.items():
        for cid in cids:
            cid_to_other_drugs[cid].add(other_drug)
    reasonable_other_drugs = set()
    for cid, other_drugs in cid_to_other_drugs.items():
        if cid in reasonable_cids:
            reasonable_other_drugs.update(other_drugs)
    return reasonable_other_drugs


def subset_to_reasonable_drugs(ranked_drugs, drug2cids, drug_annotation_dir, annot=False, out_prefix=None, out_dir=None):
    robfile = os.path.join(drug_annotation_dir, 'drug_source_moa_annotations.xlsx')
    robdf = pd.read_excel(robfile, index_col=0, dtype=str)
    reasonable_drugs = load_reasonable_drugs(drug_annotation_dir)

    rdrug_to_cids = {}
    hascid = robdf.cids.dropna()
    for drug in hascid.index:
        if drug in reasonable_drugs:
            cids = robdf.loc[drug, 'cids'].split('|')
            if cids != ['nan']:
                rdrug_to_cids[drug] = list(map(int, cids))

    reasonable_other_drugs = select_reasonable_drugs(drug2cids, rdrug_to_cids)

    cid_to_rdrug = {}
    for rdrug, cids in rdrug_to_cids.items():
        for cid in cids:
            cid_to_rdrug[cid] = rdrug

    #print([d for d in reasonable_other_drugs if d not in ranked_drugs.columns])
    reasonable_results = ranked_drugs.loc[:, reasonable_other_drugs].T

    rdrug_mechs = pd.read_excel(robfile, index_col=0, sheetname='NCI +', dtype=str).iloc[:, 7].dropna()
    rdrug_mechs.index = [str(idx).lower() for idx in rdrug_mechs.index]

    rr_mechs = pd.Series(index=reasonable_results.index)
    for drug in rr_mechs.index:
        cids = drug2cids[drug]
        for cid in cids:
            if cid in cid_to_rdrug:
                rdrug = cid_to_rdrug[cid]
                if rdrug in rdrug_mechs:
                    rr_mechs[drug] = rdrug_mechs.loc[rdrug]
                else:
                    rr_mechs[drug] = ''
    if annot:
        reasonable_results['moa'] = rr_mechs

    if out_dir is not None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        for sample in reasonable_results.columns:
            out_file = os.path.join(out_dir, '{}.{}.reasonable.annotated.csv'.format(sample, out_prefix))
            rrs = reasonable_results.loc[:, sample]
            results_plus_mech = pd.concat([rrs, rr_mechs], axis=1).sort_values(by=sample, ascending=False)
            results_plus_mech.columns = 'score moa'.split()
            results_plus_mech.index.name = 'drug'
            results_plus_mech.to_csv(out_file, float_format='%.3f')
    return reasonable_results