# INCLUDE TR01 PROJECT

## Dependencies

*  `miniforge` - check out [miniforge](https://github.com/conda-forge/miniforge). This repository holds the minimal installers for Conda and Mamba specific to conda-forge. Distribution >= Miniforge3-22.3.1-0.

* `snakemake` and other requirements

``` bash
mamba env create -f environment.yml
```

## Run analyses

Activate the conda environment

``` bash
mamba activate snakemake
```

Workflow can be executed by:

`-j` defines the number of cores or jobs to be used in parallel processing. For example, `make -j 1` would use 1 core.

``` bash
snakemake --use-conda --snakefile include_tr01.smk -j 1
```

## Datasets

### [GTEx](https://gtexportal.org/home/)

Within this repository, we obtain the GTEx Analysis V8 Gene TPMs data through the [`GTEx` Portal](https://gtexportal.org/home/), further process it, and apply [`PLIER`](https://github.com/wgmao/PLIER).
