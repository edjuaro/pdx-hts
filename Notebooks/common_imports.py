import gp
import genepattern
import os
import sys
import pandas as pd
import re
import modules.local_utils as local_utils

BASE_DIR = os.getcwd()
sys.path.append(BASE_DIR)
sys.path.append(os.path.join(BASE_DIR, 'modules'))
DATA_DIR = os.path.join(BASE_DIR, 'data')
RAW_DATA_DIR = os.path.join(DATA_DIR, 'raw')
RAW_EXP_DIR = os.path.join(RAW_DATA_DIR, 'exp')
pdx_affy_cel_dir = os.path.join(RAW_EXP_DIR, 'sebastian_3.28.17_cels/')
pdx_affy_probe_exp_file = os.path.join(RAW_EXP_DIR, 'pdx_affy_probe_exp.txt')
PDX_AFFY_ANNOT_FILE = os.path.join(RAW_EXP_DIR, '170328_PDXScreen_Anno_v5.xlsx')
GENE_ID_CONVERSION_FILE = os.path.join(RAW_DATA_DIR, 'gene_id_conversion.txt')

#preprocessed files
PREPROC_DATA_DIR = os.path.join(DATA_DIR, 'preprocessed')
PREPROC_EXP_DIR = os.path.join(PREPROC_DATA_DIR, 'exp')
PREPROC_PDX_AFFY_EXP_FILE = os.path.join(PREPROC_EXP_DIR, 'pdx_affy_exp.csv')
PREPROC_PRIMARY_AFFY_EXP_FILE = os.path.join(PREPROC_EXP_DIR, 'primary_affy_exp.csv')

#For ssGSEA
GENESETS_DIR = os.path.join(RAW_DATA_DIR, 'gene_sets')
PDX_SSGSEA_FILE = os.path.join(PREPROC_EXP_DIR, 'pdx_ssgsea.csv')
PRIMARY_SSGSEA_FILE = os.path.join(PREPROC_EXP_DIR, 'primary_ssgsea.csv')

#Used in "Process Illumina Beadchip data for 3 PDXs and primaries"
PDX_PRIMARY_BEADCHIP_SERIES_MATRIX_FILE = os.path.join(RAW_EXP_DIR, 'GSE28192_series_matrix.txt')
PDX_PRIMARY_BEADCHIP_ILMN_ANNOT_FILE = os.path.join(RAW_EXP_DIR, 'GPL6102-11574.txt') # downloaded from GEO's page for the platform
PREPROC_PDX_BEADCHIP_EXP_FILE = os.path.join(PREPROC_EXP_DIR, 'pdx_beadchip_exp.csv')
PREPROC_PRIMARY_BEADCHIP_EXP_FILE = os.path.join(PREPROC_EXP_DIR, 'primary_beadchip_exp.csv')
