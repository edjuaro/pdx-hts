
#### Before imports
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


### Before CMap
from drug_suggestion.expression.cmap import make_cmap_genesets, write_cmap_genesets
from drug_suggestion.expression.cmap import read_cmap_gct, load_cmap_drug_to_cids
from drug_suggestion.drug_annotation import subset_to_reasonable_drugs
from drug_suggestion.expression.discover import discover_from_expression, plot_discover_from_expression
from drug_suggestion.expression.controls import load_control_exp
