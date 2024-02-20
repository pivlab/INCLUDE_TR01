# plier_recount3

## Dependencies

* `conda` - for a minimal install of conda, check out [miniconda](https://docs.conda.io/en/latest/miniconda.html)	
	
* `mamba` - [drop in replacement](https://github.com/mamba-org/mamba) for conda's package manager.  Used to install `snakemake`
```
$ conda install -n base -c conda-forge mamba
```

* `snakemake` 
```
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

