# Abins mode intensities

A workflow to compute mode-resolved contributions to simulated inelastic neutron-scattering spectra.

## Introduction
This workflow calls functions included in [Mantid](https://www.mantidproject.org/), that form part of the [Abins algorithm](https://docs.mantidproject.org/nightly/algorithms/Abins-v1.html). As these functions are an "implementation detail" subject to change, the workflow is intended to execute in a Conda environment that pins Mantid to version 6.10.

The steps of the workflow are:
- Call functions within Abins to produce intensities resolved by atom and mode and write these to an sqlite database along with nuclear cross-sections from Mantid
- Perform "normal" Abins runs to produce summed fundamental and multi-phonon spectra, writing the results to .csv files
- Query the sqlite database file to collect mode-resolved sums over atoms, weighted by nuclear cross-sections, writing the results to .csv files
- Plot from .csv files to verify that results are consistent with "normal" Abins results

## Limitations
Currently the workflow will only use a single q-point of the input phonon data.

## Installation
To use this workflow, you need access to Conda (or Mamba) and Snakemake. For example, [IDAaaS](https://isis.analysis.stfc.ac.uk/) users should be able to create a suitable environment from a terminal with

```
mamba create -n snakemake -c bioconda -c conda-forge snakemake python=3.12
```

and activate it with

```
mamba activate snakemake
```


## Running the workflow

Input files are configured in *config.yml*. This workflow includes sample data for an ethanol molecule computed with GAUSSIAN16;
to use your own data and modify other workflow parameters, either edit the *config.yml* and run

```
snakemake -c 1 --sdm conda
```

or create a modified copy of *config.yml* and point to it with

```
snakemake -c 1 --sdm conda --config /path/to/my/config.yml
```
