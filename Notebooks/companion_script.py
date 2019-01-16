
#### Before imports
import gp
import genepattern
import os
import sys
import subprocess
import warnings
import re
sys.path.append('/build/')

import pandas as pd
import numpy as np
import rpy2
from sklearn.ensemble import RandomForestClassifier
from tumor_classification.medulloblastoma import classify_cavalli, classify_cho, classify_northcott
from slides import make_medullo_classification_slide, make_discover_workflow_slide, make_exp_drug_ranking_results_slide, make_intersection_slide
from utils import source_and_update_env_vars

import pickle

class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

#### Before Downloads

def log(message,extra='==>'):
    print(extra,message)
    return


def preprocess_rna_seq(setup):
    if not os.path.exists(setup.transcriptome_index_path):
        kallisto_idx_command = '{} index -i {} {}/Homo_sapiens.GRCh38.*'.format(kallisto_path, transcriptome_index_path, kallisto_dir)
        subprocess.check_call(kallisto_idx_command, shell=True).decode('utf-8')

    fq_subdirs = [os.path.join(setup.local_fastqs_dir, subdir) for subdir in os.listdir(setup.local_fastqs_dir)]
    ordered_fastqs = []
    for fq_subdir in fq_subdirs:
        if not fq_subdir.startswith('.'):
            fq_files = sorted([os.path.join(fq_subdir, file) for file in os.listdir(fq_subdir) if file.endswith('fastq.gz')])
            ordered_fastqs.extend(fq_files)
    fastqs_str = '" "'.join(ordered_fastqs)
    fastqs_str = '"'+fastqs_str+'"'

    command = '"{}" quant -i "{}" -o "{}" --bias -b 2 {}'.format(setup.kallisto_path, setup.transcriptome_index_path, setup.out_dir, fastqs_str)
    #print(command) # if you want to execute it outside this notebook in a console

    subprocess.check_call(command, shell=True)#.decode('utf-8')
    return


def run_R(command, rcript='temp_rscript.R', rlog='temp.log'):
    print("Running an R command, this may take a while and will only print output at the end of the run.")
    with open(rcript, 'w') as f:
        f.write(command)

    subprocess.call(f'Rscript {rcript} > {rlog} 2>&1 || cat {rlog}', shell=True)

    with open(rlog, 'r') as f:
        for line in f.readlines():
            print(line)
    return

def run_sleuth(setup):
    warnings.showwarning = rmagic_warning # to only print the warning text, not text + returned warning object
    from rpy2.robjects import numpy2ri
    numpy2ri.activate()
    r_command = f"""
        pkg.is.installed <- function(mypkg)
        {{
            return(mypkg %in% rownames(installed.packages()))
        }}
        source("http://bioconductor.org/biocLite.R")
        if(!pkg.is.installed('biomaRt'))
        {{
            biocLite('biomaRt')
        }}
        library("biomaRt")
        mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
          dataset = "hsapiens_gene_ensembl",
          host = 'www.ensembl.org')
        t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
        t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ext_gene = external_gene_name)

        if(!pkg.is.installed('rhdf5'))
        {{
            biocLite('rhdf5')
        }}
        if(!pkg.is.installed('devtools'))
        {{
            install.packages('devtools')
        }}
        if(!pkg.is.installed('sleuth'))
        {{
            devtools::install_github("pachterlab/sleuth")
        }}
        library(ggplot2) # req'd by sleuth
        library(dplyr, warn.conflicts = FALSE) # req'd by sleuth; masks some functions, but this doesn't affect us
        library("sleuth")
        s2c = do.call(rbind.data.frame, list(c('{setup.case_id}', '{setup.case_id}', '{setup.out_dir}')))
        colnames(s2c) <- c('sample', 'condition', 'path')
        s2c[] <- lapply(s2c, as.character)
        so <- sleuth_prep(s2c, ~condition, target_mapping=t2g, aggregation_column='ext_gene')
        write.csv(so$obs_norm[, c(2,4)], '{setup.patient_gexp_file}', quote=FALSE, row.names=FALSE)

        """
    run_R(r_command)
    numpy2ri.deactivate()
    # Return to default warning handling.
    warnings.showwarning = default_showwarning
    return

##### Before DiSCoVER
def split_discover_dataframe(df, min_score=0):
    # This is super inefficient for larger DataFrames, but let's worry about efficiency later

    out = pd.DataFrame(columns=['moa','GDSC','CTRP','CCLE','drug'])

    for index, row in df.iterrows():
        dataset = row['drug'].split('_')[0].lower()
        drug = row['drug'].split('_')[1].lower()

        if dataset == 'gdsc':
            out.loc[drug,'GDSC'] = row['score']
        elif dataset == 'ctrp':
            out.loc[drug,'CTRP'] = row['score']
        elif dataset == 'ccle':
            out.loc[drug,'CCLE'] = row['score']
        out.loc[drug,'drug'] = row['drug']
        out.loc[drug,'moa'] = row['moa']

    #     df['database'] = [drug.split('_')[0].upper() for drug in df['drug']]
    #     df['drug'] = [drug.split('_')[1].lower() for drug in df['drug']]
    #     df = df[df['score']>0]  # It feels inneficient to do this last, but Pandas complaints if I do it earlier.
    #     out = df[~df['drug'].duplicated(keep='first')].drop('score', axis=1,inplace=False)
    #     out['Drug'] = np.unique([drug.split('_')[1].lower() for drug in df['drug']])
    return out


sign_to_letter = {
    1:"+",
    -1:"-",
    '1.0':"+",
    '-1.0':"-",
    'nan':'.',
    '0.0':'.',
    '0':'.',
    0:'.',
}

def supporting_evidence(row):
    # assuming only three columns
    return sign_to_letter[str(np.sign(row[0]))]+sign_to_letter[str(np.sign(row[1]))]+sign_to_letter[str(np.sign(row[2]))]


def rank_drugs_discover(df):
    df['score'] = df.drop(['moa','drug'],axis=1,inplace=False).mean(axis=1,skipna=True).round(3)
    df['evidence'] = df.drop(['moa','score'],axis=1).apply(supporting_evidence, axis=1)
    return df.sort_values(by=['score'],ascending=False,axis=0)

###===========================================================================
### for DiSCoVER, added on 2019-01-16=========================================
###===========================================================================
def standarize_string(what):
    ascii_8bit = re.sub(r'[^ -~]', '', what).lower() # keeping only 8-bit ascii, making lowercase (http://www.catonmat.net/blog/my-favorite-regex/)
    return re.sub(r'[ _-]', '', ascii_8bit) # remove space, underscore, and dash (" _-")


# the name of the enrichment is the same as the patient ID which is stored in setup.case_id
def add_rank(df,by,name):
    df = df.sort_values(ascending=False, by=by)
    df[name] = df[by].rank(ascending=False)
    return df


# Take each disease and assign it a rank
def rank_diseases(df):
    rank_name = df.columns[-1]
    new_df = pd.DataFrame()
    to_add = pd.Series()
    new_ix = 0
    for ix,row in df.loc[:,['Disease',rank_name]].iterrows():
        for disease in row['Disease'].split('__&&__'):
            to_add['Disease'] = disease
            to_add[rank_name] = row[rank_name]
            to_add.name = new_ix
            new_df = new_df.append(to_add)
#             print(to_add)
            new_ix += 1
    return new_df


def average_disease_rank(df,rank_name):
    new_df = pd.DataFrame()
    for disease in np.unique(df['Disease']):
        new_df.loc[disease.lower(),'mean_rank'] = df[df['Disease']==disease].loc[:,rank_name].median()
        new_df.loc[disease.lower(),'weight'] = len(df[df['Disease']==disease].loc[:,rank_name])
    return new_df.sort_values(by='mean_rank')


def rank_cell_lines(setup):
    # first load some dictionaries
    cid2dic = pickle.load(file=open('discover_temp/cellosaurus/cellosaurus_cosmic_id_dic.p','rb'))
    name2dic = pickle.load(file=open('discover_temp/cellosaurus/cellosaurus_name_dic.p','rb'))

    # CCLE
    ccle = pd.read_csv(os.path.join(setup.discover_out_dir,"cell_lines_IDs_and_types_ccle.csv"),index_col=0)
    missing = 0
    for ix, row in ccle.iterrows():
        name = standarize_string(ix)
        try:
            dic = name2dic[name]
        except KeyError:
            missing += 1
            dic = {}
            dic['Names'] = ix
            dic['Disease'] = 'N/A'
            dic['CellosaurusID'] = 'N/A'
        ccle.loc[ix,'Names'] = dic['Names']
        ccle.loc[ix,'CellosaurusID'] = dic['CellosaurusID']
        ccle.at[ix,'Disease'] = dic['Disease']
    if missing >0:
        log(f'CCLE is missing {missing} cell lines (out of {len(ccle)})')
    ccle.to_excel(os.path.join(setup.discover_out_dir,"processed_cell_lines_info_ccle.xlsx"))
    # ccle.head()

    #CTRP
    ctrp = pd.read_csv(os.path.join(setup.discover_out_dir,"cell_lines_IDs_and_types_ctrp.csv"),index_col=0)
    missing = 0
    for ix, row in ctrp.iterrows():
        name = standarize_string(ix)
        try:
            dic = name2dic[name]
        except KeyError:
            missing += 1
            dic = {}
            dic['Names'] = ix
            dic['Disease'] = 'N/A'
            dic['CellosaurusID'] = 'N/A'
        ctrp.loc[ix,'Names'] = dic['Names']
        ctrp.loc[ix,'CellosaurusID'] = dic['CellosaurusID']
        ctrp.at[ix,'Disease'] = dic['Disease']
    if missing >0:
        log(f'CTRP is missing {missing} cell lines (out of {len(ctrp)})')
    ctrp.to_excel(os.path.join(setup.discover_out_dir,"processed_cell_lines_info_ctrp.xlsx"))
    # ctrp.head()

    # GDSC
    gdsc = pd.read_csv(os.path.join(setup.discover_out_dir,"cell_lines_IDs_and_types_COSMIC_IDS_gdsc.csv"),index_col=0)
    missing = 0
    for ix, row in gdsc.iterrows():
        cid = str(int(row['COSMIC ID']))
        try:
            dic = cid2dic[cid]
            gdsc.loc[ix,'Names'] = dic['Names']
            gdsc.at[ix,'Disease'] = dic['Disease']
            gdsc.loc[ix,'CellosaurusID'] = dic['CellosaurusID']
        except KeyError:
            missing+=1
    if missing>0:
        log(f"GDSC is missing {missing} cell lines (out of {len(gdsc)}). That's a bit troubling")
    gdsc.to_excel(os.path.join(setup.discover_out_dir,"processed_cell_lines_info_gdsc.xlsx"))
    # gdsc.head()

    # add ranks
    ccle = add_rank(df=ccle,by=setup.case_id,name='CCLE_rank')
    ctrp = add_rank(df=ctrp,by=setup.case_id,name='CTRP_rank')
    gdsc = add_rank(df=gdsc,by=setup.case_id,name='GDSC_rank')

    # Unique disease names, ranked
    ccle_rank = rank_diseases(df=ccle)
    ctrp_rank = rank_diseases(df=ctrp)
    gdsc_rank = rank_diseases(df=gdsc)

    #Merge the three databases
    merged = ccle_rank.rename({"CCLE_rank": "rank"}, axis='columns').append(ctrp_rank.rename({"CTRP_rank": "rank"}, axis='columns')).append(gdsc_rank.rename({"GDSC_rank": "rank"}, axis='columns'))
    merged['Disease'] = merged.apply(lambda x: x['Disease'].lower(), axis=1)

    # Only the merged rank will be printed for now
    # ccle_average_rank = average_disease_rank(df=ccle_rank,rank_name='CCLE_rank')
    # ctrp_average_rank = average_disease_rank(df=ctrp_rank,rank_name='CTRP_rank')
    # gdsc_average_rank = average_disease_rank(df=gdsc_rank,rank_name='GDSC_rank')
    merged_rank = average_disease_rank(df=merged,rank_name='rank')
    merged_rank.to_csv(os.path.join(setup.discover_out_dir,"cell_lines_rank.csv"))
    return merged_rank


###===========================================================================
### End: for DiSCoVER, added on 2019-01-16====================================
###===========================================================================


### Before CMap
from drug_suggestion.expression.cmap import make_cmap_genesets, write_cmap_genesets
from drug_suggestion.expression.cmap import read_cmap_gct, load_cmap_drug_to_cids
from drug_suggestion.drug_annotation import subset_to_reasonable_drugs
from drug_suggestion.expression.discover import discover_from_expression, plot_discover_from_expression
from drug_suggestion.expression.controls import load_control_exp


#### Before merge
def add_cmap_to_split_df(discover,cmap):
    df = discover.rename(index=str, columns={"score": "DiSCoVER", "moa":"MoA"},inplace=False)
    for index, row in cmap.iterrows():
        drug = row['drug'].lower()
        try:
            if str(df.loc[drug,'CMAP']) != 'nan':
                #drug already exists and a CMAP score has been added
                df.loc[drug,'CMAP'] = np.nanmean([row['score'], df.loc[drug,'CMAP']])
            else:
                # drug already exists but a CMAP score has not been added
                df.loc[drug,'CMAP'] = row['score']
        except KeyError:
            # drug didn't exist therefore a CMAP score had not been added
            df.loc[drug,'CMAP'] = row['score']
            df.loc[drug,'MoA'] = row['moa']
        if str(df.loc[drug,'evidence'])=='nan': #numpy is not letting me use np.isnan()
            df.loc[drug,'evidence'] = '...'+sign_to_letter[str(np.sign(row['score']))]
        else:
            df.loc[drug,'evidence'] = str(df.loc[drug,'evidence'])+sign_to_letter[str(np.sign(row['score']))]

    #update those rows which are not in cmap
    for index, row in df.iterrows():
        if len(row['evidence'])==3:
            df.loc[index,'evidence']  = df.loc[index,'evidence']+'.'

    return df[['drug','MoA','GDSC','CTRP','CCLE','DiSCoVER','CMAP','evidence']]

def rank_combined_df(df):
#     df['average'] = df.drop(['MoA','GDSC','CTRP','CCLE','support'],axis=1,inplace=False).mean(axis=1,skipna=True).round(3)
    df.sort_values(by=['DiSCoVER'],ascending=False,axis=0,inplace=True)
    df['DiSCoVER rank'] = range(1, len(df) + 1)
    df.sort_values(by=['CMAP'],ascending=False,axis=0,inplace=True)
    df['CMAP rank'] = range(1, len(df) + 1)
    df['Average rank'] = df[['DiSCoVER rank','CMAP rank']].mean(axis=1)
    df.sort_values(by=['Average rank'],ascending=True,axis=0,inplace=True)
    df = df[((df['DiSCoVER']>0.001).values) & ((df['CMAP']>0.001).values)]
#     combined_df = combined_df[((combined_df['DiSCoVER']>0.001).values) & ((combined_df['CMAP']>0.001).values)]

    return df.drop(['GDSC','CTRP','CCLE'],axis=1,inplace=False)

## For discover
from collections import defaultdict
import utils
drug_annotation_dir = os.path.dirname('/build/drug_suggestion/drug_annotation/')
sys.path.append(os.path.join(drug_annotation_dir))


def load_reasonable_drugs():
    reasonable_drugs_file = os.path.join(drug_annotation_dir, 'clinically_relevant_drugs.csv')
    with open(reasonable_drugs_file, 'r') as f:
        rdrugs = [row.strip() for row in f.readlines()]
    return rdrugs


reasonable_drugs = load_reasonable_drugs()


def select_reasonable_drugs(other_drug_to_cids, rdrug_to_cids):
    reasonable_cids = utils.all_unique_dict_values(rdrug_to_cids)
    cid_to_other_drugs = defaultdict(set)
    for other_drug, cids in other_drug_to_cids.items():
        for cid in cids:
            cid_to_other_drugs[cid].add(other_drug)
    reasonable_other_drugs = set()
    for cid, other_drugs in cid_to_other_drugs.items():
        if cid in reasonable_cids:
            reasonable_other_drugs.update(other_drugs)
    return reasonable_other_drugs


def format_drugs(ranked_drugs, drug2cids, annot=False, out_prefix=None, out_dir=None):
    robfile = os.path.join(drug_annotation_dir, 'drug_source_moa_annotations.xlsx')
    robdf = pd.read_excel(robfile, index_col=0, dtype=str)

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

    # The line below was commented out on 2018-12-14
    # reasonable_results = ranked_drugs.loc[:, reasonable_other_drugs].T
    reasonable_results = ranked_drugs.T

    rdrug_mechs = pd.read_excel(robfile, index_col=0, sheetname='NCI +', dtype=str).iloc[:, 7].dropna()
    rdrug_mechs.index = [str(idx).lower() for idx in rdrug_mechs.index]

    rr_mechs = pd.Series(index=reasonable_results.index)
    # counter1 = 0
    # counter2 = 0
    for drug in rr_mechs.index:
        try:
            cids = drug2cids[drug]
            # counter1 += 1
        except KeyError:
            cids = 'Not Clinically Relevant'
            # counter2 += 1
        for cid in cids:
            if cid in cid_to_rdrug:
                rdrug = cid_to_rdrug[cid]
                if rdrug in rdrug_mechs:
                    rr_mechs[drug] = rdrug_mechs.loc[rdrug]
                else:
                    rr_mechs[drug] = ''
            else:
                rr_mechs[drug] = 'Not Clinically Relevant'
    # print(counter1,counter2)
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
    return results_plus_mech
