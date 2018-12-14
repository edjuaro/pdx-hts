import os
import json
import sys

import pandas as pd
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
from joblib import Parallel, delayed

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from modules.controls import make_dx_disease_gdict
from modules.gsea import ssgsea
from modules.information import compute_ic

supported_cls = 'gdsc ctrp ccle'.split()


def load_discover_drug_to_cids(discover_data_dir, cls):
    d = {}
    for cl in cls:
        if cl not in supported_cls:
            raise ValueError('{} not a supported cl; try one of {}'.format(cl, supported_cls))
        d.update(load_cl_drug_to_cids(discover_data_dir, cl, prefix=True))
    return d


def load_cl_drug_to_cids(discover_data_dir, cl, prefix=False):
    to_pubchem_cid_file = os.path.join(discover_data_dir, cl, '{}_name_to_pubchem_cid.json'.format(cl))
    to_cid = json.load(open(to_pubchem_cid_file, 'r'))
    if prefix:
        return {'{}_{}'.format(cl, key): value for key, value in to_cid.items()}
    return to_cid


def make_discover_genesets(exp, control_exp, cl_exp=None, method='rankdif'):
    if cl_exp is None:
        valid_genes = None
    else:
        valid_genes = cl_exp.columns.tolist()
    gsets = make_dx_disease_gdict(exp, control_exp, n_genes=150, up=True, method=method, valid_genes=valid_genes, sets=False)
    return gsets


def discover_from_signature(discover_data_dir, disease_gdict, alpha=0.75, verbose=False):
    return discover(discover_data_dir, disease_gdict=disease_gdict, alpha=alpha, verbose=verbose)


def discover_from_expression(discover_data_dir, exp, control_exp, alpha=0.75, verbose=False):
    return discover(discover_data_dir, exp=exp, control_exp=control_exp, alpha=alpha, verbose=verbose)


def discover(discover_data_dir, exp=None, control_exp=None, disease_gdict=None, alpha=0.75, verbose=False):
    cl_names = 'ctrp gdsc ccle'.split()
    all_ics = []
    for cl_name in cl_names:
        if verbose:
            print("Loading {}".format(cl_name.upper()))
        hdf_file = os.path.join(discover_data_dir, cl_name, 'store.h5')
        cl_exp = pd.read_hdf(hdf_file, 'exp')
        cl_dr = pd.read_hdf(hdf_file, 'dr')
        # they are already reduced to cell lines they have in common
        # exp_df, dr_df = reduce_to_common_idxs([exp_df, dr_df], axis=1)
        if cl_name in ['gdsc', 'ctrp']:
            cl_dr *= -1
        if verbose:
            print('Projecting disease gene sets')
        if disease_gdict is None:
            disease_gdict = make_discover_genesets(exp, control_exp, cl_exp=cl_exp)
        ssgsea_df = ssgsea(cl_exp, disease_gdict, alpha=alpha, n_jobs=-1)
        if verbose:
            print('Matching to drug response profiles')
        ics = compute_discover_ics(ssgsea_df, cl_dr, n_jobs=-1)
        ics.columns = ['{}_{}'.format(cl_name, idx) for idx in ics.columns]
        all_ics.append(ics)
    combined_results = pd.concat(all_ics, axis=1)
    return combined_results


def compute_discover_ics(ssgsea_df, dr_df, n_jobs=1):
    if n_jobs == -1:
        n_jobs = os.cpu_count()
    diseases = ssgsea_df.columns
    drugs = dr_df.columns
    n_diseases = len(diseases)
    n_drugs = len(drugs)
    ic_inputs = [(ssgsea_df.loc[:, disease], dr_df.loc[:, drug]) for disease in diseases for drug in drugs]
    ics_list = Parallel(n_jobs=n_jobs)(delayed(compute_ic)(*inputs) for inputs in ic_inputs)
    ics = pd.DataFrame(np.array(ics_list).reshape((n_diseases, n_drugs)), index=diseases, columns=drugs)
    return ics


def plot_discover_from_signature(sample_name, discover_results, disease_gdict, cl='ctrp', alpha=0.75, out_file=None, min_nonnull_frac=0.5):
    return plot_discover(sample_name, discover_results, disease_gdict=disease_gdict, cl=cl, alpha=alpha, out_file=out_file, min_nonnull_frac=min_nonnull_frac)


def plot_discover_from_expression(sample_name, discover_results, exp, control_exp, cl='ctrp', alpha=0.75, out_file=None, min_nonnull_frac=0.5):
    return plot_discover(sample_name, discover_results, exp=exp, control_exp=control_exp, cl=cl, alpha=alpha, out_file=out_file, min_nonnull_frac=min_nonnull_frac)


def plot_discover(discover_data_dir, sample_name, discover_results, disease_gdict=None, exp=None, control_exp=None, cl='ctrp', alpha=0.75, out_file=None, min_nonnull_frac=0.5):
    hdf_file = os.path.join(discover_data_dir, cl, 'store.h5')
    cl_exp = pd.read_hdf(hdf_file, 'exp')
    cl_dr = pd.read_hdf(hdf_file, 'dr')
    if cl in ['ctrp', 'gdsc']:
        cl *= -1

    if disease_gdict is None:
        disease_gdict = make_discover_genesets(exp, control_exp, cl_exp=cl_exp)

    sub_disease_gict = {sample_name: disease_gdict[sample_name]}
    signature = ssgsea(cl_exp, sub_disease_gict, alpha=alpha).loc[:, sample_name]

    dr_ma = ma.array(cl_dr.values, mask=cl_dr.isnull())
    ma_mean = dr_ma.mean(axis=0)
    ma_std = dr_ma.std(axis=0)
    ma_scaled = ma.divide(ma.subtract(dr_ma, ma_mean), ma_std)

    normed_dr = pd.DataFrame(ma_scaled, index=cl_dr.index, columns=cl_dr.columns)

    drug_nonnull_frac = normed_dr.notnull().sum(axis=0) / normed_dr.shape[0]
    enough_nonnull = drug_nonnull_frac > min_nonnull_frac
    viab_df = normed_dr.T[enough_nonnull].T

    ctrp_disco = discover_results.loc[:, [d for d in discover_results.columns if d.startswith('ctrp_')]]

    ctrp_disco.columns = [col.split('_')[1] for col in ctrp_disco.columns]
    disco_passed = ctrp_disco.loc[:, enough_nonnull[enough_nonnull].index]

    ranked_drugs = disco_passed.loc[sample_name].sort_values(ascending=False)
    sns.set(font_scale=2.0)
    fig = plt.figure()
    fig.set_size_inches(7.5, 10)
    gs = gridspec.GridSpec(10, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1:, 0])
    sorted_sig = signature.sort_values(ascending=False)
    n_cell_lines = len(sorted_sig)

    ndrugs = 10
    top_drugs = list(ranked_drugs.head(ndrugs).index)
    bottom_drugs = list(ranked_drugs.tail(ndrugs).index)
    featured_drugs = top_drugs + bottom_drugs

    xr = range(n_cell_lines)
    y = sorted_sig
    ax1.plot(xr, y)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_xlim(left=0, right=n_cell_lines)
    ax1.grid(False)
    ax1.set_ylabel('{}\nsignature\nenrichment'.format(sample_name), rotation=0, fontsize=14)
    ax1.yaxis.set_label_coords(-0.125, 0.05)
    ax1.yaxis.set_label_position("left")
    # ax1.set_facecolor('white')
    cbar_ax = fig.add_axes([.92, 0.125, .03, 0.68])
    toplot = viab_df.loc[sorted_sig.index, featured_drugs].T  # .clip(lower=-2, upper=2)
    toplot = toplot.rename(index={idx: idx.replace('hydrochloride', '') for idx in toplot.index})
    hm = sns.heatmap(toplot, ax=ax2, cbar_ax=cbar_ax, cbar_kws={'label': 'normalized viability'})
    yt = cbar_ax.get_yticks()
    # cbar_ax.set_yticks([yt[0], 1, 2, 3, 4])
    cbar_ax.set_yticks([-1, 0, 1, 2, 3])
    cbar_ax.set_yticklabels(['low', '', '', '', 'high'])
    hm.set_xticklabels([])
    hm.set_xlabel('Ranked cell lines')
    hm.set_ylabel('Ranked drugs (top and bottom {})'.format(ndrugs))
    for i, tl in enumerate(hm.get_yticklabels()):
        tl.set_color('b' if i < ndrugs else 'r');  # since yticklabels are from bottom up
    if out_file is not None:
        plt.savefig(out_file, bbox_inches='tight')