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
        f'{config.output}/gtex/GTEx_v8_gene_median_tpm.rds',
        f'{config.output}/gtex/gtex_plier.rds',
        expand(f'{config.output}/gtex/robustness/gtex_plier_rob_{{replicate}}.rds', replicate=range(1, 11))


# INSTALL DEPENDENCIES
rule install_PLIER:
    """
    Download and install PLIER inside the conda env
    """
    output:
        f'{config.logs}/PLIER_installed.txt'
    conda:
        'envs/gtex.yaml'
    shell:
        """
        R -e "devtools::install_github('wgmao/PLIER')" && touch {output}
        """

# DOWNLOAD DBs
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

# PROCESS AND ANALIZE DATA
rule process_gtex_data:
    """
    Processes GTEx gene expression TSV file, renames columns, and saves it in RDS format
    """
    input:
        script = "scripts/process_gtex_data.R",
        gtex_data = f'{rules.download_gtex_data.output.gtex_data}'
    output:
        gtex_data_p = f'{config.output}/gtex/GTEx_v8_gene_median_tpm.rds', 
    conda:
        'envs/gtex.yaml',
    shell:
        """
        Rscript {input.script} {input.gtex_data} {output.gtex_data_p}
        """


rule plier_gtex:
    """
    Processes GTEx gene expression data using PLIER for pathway analysis, including data Z-score normalization
    """
    input:
        script = "scripts/plier_gtex.R",
        gtex_data_p = f'{rules.process_gtex_data.output.gtex_data_p}'
    output:
        gtex_plier = f'{config.output}/gtex/gtex_plier.rds', 
    conda:
        'envs/gtex.yaml',
    shell:
        """
        Rscript {input.script} {input.gtex_data_p} {output.gtex_plier}
        """

# rule analyze_plier_gtex:
#     """
#     """
#     input:
#         script = "scripts/render_nbs.sh",
#         input_nb = "nbs/10_gtex/GTEx_PLIER_exploration.ipynb"
#         gtex_plier = f'{config.output}/gtex/gtex_plier.rds', 
#     output:
#         _ = "", 
#     conda:
#         'envs/jupyter.yaml',
#     shell:
#         """
#         {input.script} -p INPUT_PLIER_MODEL_FILE {output.plier_gtex.gtex_plier}
#         """
    
rule plier_gtex_robustness:
    """
    Processes GTEx gene expression data using PLIER for pathway analysis, including data Z-score normalization 10 times to assay the robustness of the results
    """
    input:
        script = "scripts/plier_gtex.R",
        gtex_data_p = f'{rules.process_gtex_data.output.gtex_data_p}'
    output:
        gtex_plier = f'{config.output}/gtex/robustness/gtex_plier_rob_{{replicate}}.rds'
    conda:
        'envs/gtex.yaml'
    shell:
        """
        Rscript {input.script} {input.gtex_data_p} {output.gtex_plier}
        """
