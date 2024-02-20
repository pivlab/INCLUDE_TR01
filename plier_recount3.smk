"""
PLIER recount3
"""

# ==============================================================================
# IMPORT LIBS
# ==============================================================================

from os.path import basename
from glob import glob
from types import SimpleNamespace
import os
import utils
import pandas as pd


# ==============================================================================
# SETUP
# ==============================================================================

configfile: 'conf/config.yaml'

# This allows for dot (.) access instead of dict access
# of configuration parameters.
# EXAMPLE: config['ABC'] becomes config.ABC

config = SimpleNamespace(**config)

# ==============================================================================
# RULES
# ==============================================================================
rule All:
    input:
        # Download data
        f'{config.data}/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',
        # Install dependencies
        f'{config.logs}/PLIER_installed.txt',
        # Process files
        f'{config.output}/gtex/GTEx_v8_gene_median_feather'


# INSTALL DEPENDENCIES
rule install_PLIER:
    output:
        f'{config.logs}/PLIER_installed.txt'
    conda:
        'envs/gtex.yaml'
    shell:
        """
        R -e "devtools::install_github('wgmao/PLIER')" && touch {output}
        """

rule download_gtex_data:
    """
    Download GTEx v8 expression data
    """
    input:
        script = "scripts/download_gtex_data.sh",
    output:
        gtex_data = f'{config.data}/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct', 
    shell:
        """
        {input.script} {config.gtex_link} {output.gtex_data}
        """

rule process_gtex_data:
    """
    Read GTEx v8 data and save the entire gene expression data
    """
    input:
        script = "scripts/process_gtex_data.R",
        gtex_data = f'{rules.download_gtex_data.output.gtex_data}'
    output:
        gtex_data_p = f'{config.output}/gtex/GTEx_v8_gene_median_feather', 
    conda:
        'envs/gtex.yaml',
    shell:
        """
        Rscript {input.script} {input.gtex_data} {output.gtex_data_p}
        """
