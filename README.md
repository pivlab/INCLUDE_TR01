# plier_recount3

## Dependencies

* `conda` - for a minimal install of conda, check out [miniconda](https://docs.conda.io/en/latest/miniconda.html)	
	
* `mamba` - [drop in replacement](https://github.com/mamba-org/mamba) for conda's package manager.  Used to install `snakemake`

``` bash
conda install -n base -c conda-forge mamba
```

* `snakemake` 

``` bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Install python requirments

``` bash
conda activate snakemake
pip install -r requirements.txt
```

## Run analyses

Activate the conda environment

``` bash
conda activate snakemake
```

Workflow can be executed by:

`-j` defines the number of cores or jobs to be used in parallel processing. For example, `make -j 1` would use 1 core.

``` bash
snakemake --use-conda --snakefile plier_recount3.smk -j 1
```

