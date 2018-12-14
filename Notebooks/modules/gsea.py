import os
import itertools

import numpy as np
import pandas as pd
from joblib import Parallel, delayed


def consecutive_pairs(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def _base_gsea(ranked_genes, sets_to_genes, collect_func, alpha=0.75):
    """
    Basic idea:
    -make weighted ecdf for hits
    -make ecdf for misses
    -take elementwise difference
    -collect_func() to get result (max_abs() for gsea, np.sum() for ssgsea)

    Might still be able to speed up using array funcs rather than iterating?
    """
    n_genes = len(ranked_genes)
    ranks = list(range(n_genes))
    gene_to_rank = dict(zip(ranked_genes, ranks))
    enrichment_scores = {}
    for set_name, set_genes in sets_to_genes.items():
        ranked_set_genes = [gene for gene in set_genes if gene in gene_to_rank]
        n_non_set_genes = float(n_genes - len(ranked_set_genes))
        hit_ranks = [gene_to_rank[gene] for gene in ranked_set_genes]
        misses = np.ones_like(ranks)
        misses[hit_ranks] = 0
        cum_misses = np.cumsum(misses)
        miss_ecdf = cum_misses / n_non_set_genes
        cum_hits = np.zeros_like(ranks)
        if len(hit_ranks) > 0:
            cum_hit_sum = 0
            sorted_hit_ranks = sorted(hit_ranks)
            # add one so ranks to weight start from 1, not zero
            # however, it's convenient to start at zero otherwise so I can index using the ranks
            weighted_ranks = (np.array(sorted_hit_ranks) + 1) ** alpha
            hit_rank_pairs = list(consecutive_pairs(sorted_hit_ranks))  # given [a, b, c, d] yields (a, b), (b, c), (c, d)
            for i, (idx1, idx2) in enumerate(hit_rank_pairs):
                cum_hit_sum += weighted_ranks[i]
                cum_hits[idx1:idx2] = cum_hit_sum
            cum_hit_sum += weighted_ranks[-1]
            cum_hits[sorted_hit_ranks[-1]:] = cum_hit_sum
            weighted_hit_ecdf = cum_hits / cum_hit_sum
        else:
            weighted_hit_ecdf = cum_hits  # still np.zeros_like(ranks)
        ecdf_dif = np.subtract(weighted_hit_ecdf, miss_ecdf)
        enrichment_score = collect_func(ecdf_dif)
        enrichment_scores[set_name] = enrichment_score
    return pd.Series(enrichment_scores)


def ssgsea(exp_data, sets_to_genes, alpha=0.75, n_jobs=-1):
    """
    Single-sample GSEA as described in Barbie et al. (2009)
    :param exp_df: Pandas DataFrame or Series of expression values, (n_samples, n_genes) or (n_genes,)
    :param sets_to_genes: dictionary with set names as keys and sets of genes as values, e.g. {'set1': {'g1', 'g2'}}
    :param alpha: float, weighting factor between zero and one.
        Smaller values give more weight to top/bottom of list.
    :return: Pandas DataFrame or Series of expression projected onto gene sets
    """
    if isinstance(exp_data, pd.Series):
        return ssgsea_per_sample(exp_data, sets_to_genes, alpha=alpha)
    elif isinstance(exp_data, pd.DataFrame):
        if n_jobs == -1:
            n_jobs = os.cpu_count()
        ssgsea_inputs = [(exp_data.loc[sample], sets_to_genes, alpha) for sample in exp_data.index]
        ssgseas = Parallel(n_jobs=n_jobs)(delayed(ssgsea_per_sample)(*inputs) for inputs in ssgsea_inputs)
        return pd.concat(ssgseas, axis=1).T
        #return exp_data.apply(ssgsea_per_sample, axis=1, args=(sets_to_genes, alpha))
    else:
        raise ValueError("exp_data must be Pandas DataFrame or Series")


def ssgsea_per_sample(exp_series, sets_to_genes, alpha):
    sorted_exp_series = exp_series.sort_values(ascending=False)
    enrichment_scores = _base_gsea(sorted_exp_series.index, sets_to_genes, collect_func=np.sum, alpha=alpha)
    enrichment_scores /= 0.5 * len(exp_series)  # maximum possible score
    enrichment_scores.name = exp_series.name
    return enrichment_scores