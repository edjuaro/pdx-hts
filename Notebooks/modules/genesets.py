import os
import itertools
from collections import defaultdict


def load_genesets(genesets_dir, msigdb_version='6.1', which='all'):
    msigdb_all_genesets_gmt = os.path.join(genesets_dir, 'msigdb.v{}.symbols.gmt'.format(msigdb_version))
    ipa_gmt = os.path.join(genesets_dir, 'IPA_regulators.gmt')
    myc_gmt = os.path.join(genesets_dir, 'MYC_signatures.gmt')
    c1_gmt = os.path.join(genesets_dir, 'c1.all.v{}.symbols.gmt'.format(msigdb_version))
    c2_gmt = os.path.join(genesets_dir, 'c2.all.v{}.symbols.gmt'.format(msigdb_version))
    c3_gmt = os.path.join(genesets_dir, 'c3.all.v{}.symbols.gmt'.format(msigdb_version))
    c4_gmt = os.path.join(genesets_dir, 'c4.all.v{}.symbols.gmt'.format(msigdb_version))
    c5_gmt = os.path.join(genesets_dir, 'c5.all.v{}.symbols.gmt'.format(msigdb_version))
    bp_gmt = os.path.join(genesets_dir, 'c5.bp.v{}.symbols.gmt'.format(msigdb_version))
    c6_gmt = os.path.join(genesets_dir, 'c6.all.v{}.symbols.gmt'.format(msigdb_version))
    c7_gmt = os.path.join(genesets_dir, 'c7.all.v{}.symbols.gmt'.format(msigdb_version))
    c3_tft_gmt = os.path.join(genesets_dir, 'c3.tft.v{}.symbols.gmt'.format(msigdb_version))
    c3_mir_gmt = os.path.join(genesets_dir, 'c3.mir.v{}.symbols.gmt'.format(msigdb_version))
    all_genesets_gmt = os.path.join(genesets_dir, 'all_genesets.gmt')
    hallmark_gmt = os.path.join(genesets_dir, 'h.all.v{}.symbols.gmt'.format(msigdb_version))
    complementary_signatures_gmt = os.path.join(genesets_dir, 'complementary.signatures.v9.gmt')
    biocarta_gmt = os.path.join(genesets_dir, 'c2.cp.biocarta.v{}.symbols.gmt'.format(msigdb_version))
    kegg_gmt = os.path.join(genesets_dir, 'c2.cp.kegg.v{}.symbols.gmt'.format(msigdb_version))
    reactome_gmt = os.path.join(genesets_dir, 'c2.cp.reactome.v{}.symbols.gmt'.format(msigdb_version))
    pid_gmt = os.path.join(genesets_dir, 'c2.cp.pid.v{}.symbols.gmt'.format(msigdb_version))
    cgp_gmt = os.path.join(genesets_dir, 'c2.cgp.v{}.symbols.gmt'.format(msigdb_version))
    cp_gmt = os.path.join(genesets_dir, 'c2.cp.v{}.symbols.gmt'.format(msigdb_version))
    arg_to_gmt_dict = {    'ipa': ipa_gmt,
                           'myc': myc_gmt,
                        'msigdb': msigdb_all_genesets_gmt,
                      'hallmark': hallmark_gmt,
                      'comp_sig': complementary_signatures_gmt,
                            'c1': c1_gmt,  # positional
                            'c2': c2_gmt,  # curated pathways, includes Kegg, Reactome, BioCarta, etc
                          'kegg': kegg_gmt,
                      'biocarta': biocarta_gmt,
                      'reactome': reactome_gmt,
                           'pid': pid_gmt,
                           'cgp': cgp_gmt,
                            'cp': cp_gmt,
                            'c3': c3_gmt,  # regulatory motifs
                            'c4': c4_gmt,  # computationally-derived cancer signatures
                            'c5': c5_gmt,  # GO
                            'bp': bp_gmt,
                            'c6': c6_gmt,  # oncogenic signatures
                            'c7': c7_gmt,  # immunologic signatures
                           'tft': c3_tft_gmt,
                           'mir': c3_mir_gmt,
                           'all': all_genesets_gmt
                      }
    geneset_names_string = ", ".join(arg_to_gmt_dict.keys())
    if isinstance(which, str):
        if which not in arg_to_gmt_dict:
            raise ValueError("{} not a supported geneset; try one of {}".format(which, geneset_names_string))
        p2g = read_gmt(arg_to_gmt_dict[which])
        p2g = fix_synonyms(p2g)
    elif type(which) is list:
        all_items = [load_genesets(which=wh).items() for wh in which]
        p2g = dict(itertools.chain(*all_items))
    return p2g


def fix_synonyms(p2g, take='both'):
    new_p2g = defaultdict(set)
    for pathway, genes in p2g.items():
        for gene in genes:
            synonyms = gene.split(" /// ")
            if len(synonyms) > 1:
                if take == 'both':
                    new_p2g[pathway].update(synonyms)
                elif take == 'first':
                    new_p2g[pathway].add(synonyms[0])
            else:
                new_p2g[pathway].add(gene)
    return p2g


def write_gmt(sets_to_genes, gmt_file, append=False):
    strings = []
    for s, genes in sets_to_genes.items():
        try:
            string = '{}\t\t{}'.format(s, '\t'.join(map(str, genes)))
        except Exception as e:
            print(s)
            print(genes)
            raise(e)
        strings.append(string)
    with open(gmt_file, 'a' if append else 'w') as f:
        f.write('\n'.join(strings))
        f.write('\n')


def read_gmt(gmt_file):
    with open(gmt_file, 'r') as f:
        lines = f.readlines()
    p2g = {}
    for line in lines:
        entries = line.strip().split("\t")
        pathway_name = entries[0]
        genes = set(entries[2:])
        if len(genes) > 0:
            p2g[pathway_name] = genes
    return p2g
