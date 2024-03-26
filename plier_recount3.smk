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
        ### GTEx dataset ###
        # Download data
        #f'{config.data}/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',
        # Install dependencies
        #f'{config.logs}/PLIER_installed.txt',
        # Process files
        #f'{config.output}/gtex/GTEx_v8_gene_median_tpm.rds',
        #f'{config.output}/gtex/gtex_plier_prep.rds',
        # run GTEx with multiple parameters
        # expand(f'{config.output}/gtex/plier_result_k{{parameter_k}}_frac{{frac}}.rds',
        #       parameter_k=[0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2],
        #      frac=[0.25, 0.5, 0.7, 1]),
        #expand(f'{config.output}/gtex/robustness/gtex_plier_rob_{{replicate}}.rds', replicate=range(1, 11)),
        # Jupyter notebooks
        "nbs/10_gtex/GTEx_PLIER_exploration.run.txt",
        "nbs/10_gtex/GTEx_PLIER_robustness.run.txt",
        "docs/index.html"


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
rule gtex_data_preprocess:
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


rule gtex_plier_prep:
    """
    Processes GTEx gene expression data including Z-score normalization
    """
    input:
        script = "scripts/plier_preprocess.R",
        gtex_data_p = f'{rules.gtex_data_preprocess.output.gtex_data_p}'
    output:
        gtex_plier_prep = f'{config.output}/gtex/gtex_plier_prep.rds', 
    conda:
        'envs/gtex.yaml',
    shell:
        """
        Rscript {input.script} {input.gtex_data_p} {output.gtex_plier_prep}
        """

rule gtex_plier_run:
    """
    Processes GTEx gene expression data using PLIER with default k and frac parameters.
    """
    input:
        script = "scripts/plier_run.R",
        gtex_data_prep = f'{rules.gtex_plier_prep.output.gtex_plier_prep}'
    output:
        gtex_plier = f'{config.output}/gtex/plier_result_k{{parameter_k}}_frac{{frac}}.rds'
    params:
        parameter_k = "{parameter_k}",
        frac = "{frac}"
    conda:
        'envs/gtex.yaml',
    shell:
        """
        Rscript {input.script} {input.gtex_data_prep} {output.gtex_plier} {params.parameter_k} {params.frac}
        """

rule plier_gtex_robustness:
    """
    Processes GTEx gene expression data using PLIER for pathway analysis, including data Z-score normalization 10 times to assay the robustness of the results
    """
    input:
        script = "scripts/plier_gtex.R",
        gtex_data_p = f'{rules.gtex_plier_run.output.gtex_plier}'
    output:
        gtex_plier = f'{config.output}/gtex/robustness/gtex_plier_rob_{{replicate}}.rds'
    conda:
        'envs/gtex.yaml'
    shell:
        """
        Rscript {input.script} {input.gtex_data_p} {output.gtex_plier}
        """

rule GTEx_PLIER_exploration_analysis:
    """
    Analysis of the robustness of the processed GTEx gene expression data using PLIER for pathway analysis
    """
    input:
        script = "scripts/render_nbs.sh",
        input_nb = "nbs/10_gtex/GTEx_PLIER_exploration.ipynb",
        gtex_plier = "output/gtex/plier_result_k1_frac0.7.rds"
    output:
        "nbs/10_gtex/GTEx_PLIER_exploration.run.txt"
    conda:
        'envs/gtex.yaml',
    shell:
        """
        {input.script} {input.input_nb} -p INPUT_PLIER_MODEL_FILE {input.gtex_plier}
        """

rule GTEx_PLIER_robustness_analysis:
    """
    Analysis of the robustness of the processed GTEx gene expression data using PLIER for pathway analysis
    """
    input:
        script = "scripts/render_nbs.sh",
        input_nb = "nbs/10_gtex/GTEx_PLIER_robustness.ipynb",
        gtex_data_p = os.path.dirname(f'{rules.plier_gtex_robustness.output.gtex_plier}') 
    output:
        "nbs/10_gtex/GTEx_PLIER_robustness.run.txt"
    conda:
        'envs/gtex.yaml',
    shell:
        """
        {input.script} {input.input_nb} -p INPUT_PLIER_ROB_DIR {input.gtex_data_p}
        """

rule render_website:
    """
    Render website using Quarto
    """
    output:
        "docs/index.html"
    shell:
        """
        quarto render
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
    

# rule plier_gtex_robustness_rseed:
#     """
#     Processes GTEx gene expression data using PLIER for pathway analysis, including data Z-score normalization 10 times to assay the robustness of the results
#     modifying the rseed
#     """
#     input:
#         script = "scripts/plier_gtex.R",
#         gtex_data_p = f'{rules.process_gtex_data.output.gtex_data_p}'
#     output:
#         gtex_plier = f'{config.output}/gtex/robustness/gtex_plier_rob_{{replicate}}.rds'
#     conda:
#         'envs/gtex.yaml'
#     shell:
#         """
#         Rscript {input.script} {input.gtex_data_p} {output.gtex_plier}
#         """