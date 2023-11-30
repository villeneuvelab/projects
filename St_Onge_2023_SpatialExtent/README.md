# St-Onge et al. (2023) - Spatial Extent Project

**Author**: Frederic St-Onge
**Date**: September 7th 2023

The current folder contains all the code necessary to reproduce the 
analyses in our paper. To match the Python environment, you can recreate 
the conda environment locally from the `environment.yml` file if you have 
conda installed on your computer by following the commands:

```bash
$ conda env create -f environment.yml
$ conda activate spex_env
$ conda env list #Verifies the installation
```

Alternatively, if you do not have conda, you can simply create a virtual 
environment with `pip` and install `sihnpy`. Make sure to do this using 
Python 3.9 to match the software version. Installing `sihnpy` will also 
pull `pandas` and `numpy` which are both used to clean the data.

Most of the statistical analyses were done in the R markdown
document `src/r/spex_analyses.Rmd`.