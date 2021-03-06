# Lane, Bovy, & Mackereth (2021)

This repository contains code to replicate the results of our [paper](https://ui.adsabs.harvard.edu/abs/2021arXiv210609699L/abstract), but also to perform similar analyses. The `master` branch contains notebooks which have been cleaned,  the `run-nbs` branch contains notebooks which have been run sucessfully, to give an idea of what the results should look like. We encourage others to use this as a template to analyze their own data.

## Requirements

The normal `scipy`/`numpy` stack as well as `matplotlib` and `astropy`. We recommend the development version of [galpy](https://github.com/jobovy/galpy), as that may contain recent updates to the spherical distribution function code. Version 1.7 should suffice though. If you need to download APOGEE data we use [apogee](https://github.com/jobovy/apogee) and if you need to download Gaia data we use [gaia_tools](https://github.com/jobovy/gaia_tools). 

Notebooks 1 (`notebooks/1-get_data/get_data.ipynb`) and 2 (`notebooks/prepare_data.ipynb`) below will download and prepare the data we use in our paper. We encourage those with other data sets to look at those notebooks to get an idea of the kind of data that we use, and how we prepare it for use with the distribution functions. Notebook 3 (`notebooks/sample_dfs`) shows how the data are used in combination with the distribution functions to generate mock catalogues of 

## The data

We use [APOGEE DR16](https://www.sdss.org/dr16/) and *Gaia* DR2 data in our paper. The first notebook `notebooks/get_data.ipynb` will download these data and prepare them, as well as auxilliary data components such as the selection function and statistical sample. This will take awhile (hours) to run.

The second notebook, `notebooks/2-process_data/process_data.ipynb`, will take the raw data, extract a high-quality subsample, and divide it into subgroups of stars (thin disk/thick disk/halo) for the final analysis. Note that this is specifically the approach to matching stellar populations with distribution functions that we took. Others with different data (or even APOGEE data) should feel free to divide their samples into stellar populations as they see fit. It's just important to consider that most stellar samples will contain a mix of stars from different stellar populations, each with their own distribution function. 

## The distribution functions

The third step is in scripts, `scripts/3-sample_dfs/*` will take the high-quality data from notebook 2 (divided into thin disk/thick disk/halo) and use distribution functions to generate kinematic samples. The process is divided into multiple scripts since the process can be time consuming. The scripts can be run in parallel, they don't depend on one another.

## The figures

`notebooks/4-figures/*.ipynb` generate the different figures for the paper.
