This is the software repository for the paper "Influence of European heat waves on boundary layer heights" by Till Fohrmann, Svenja Szemkus, Arianna Valmassoi, and Petra Friederichs, which is being submitted to Journal of Geophysical Research: Atmospheres.

The contents are:
- "data": Containing all the data for evaluation, which is accessible from:
- "figs": Output location for the figures shown in the paper.
- "download_era5_blh.py": Download script for accessing the ERA5 BLH diagnostic as we used it in this study
- "emd_bootstrap.ipynb": A jupyter notebook which handels the grid point wise hypothesis testing using Julia (see .toml for environment specifications)
- "fdbck_tools.py": Some functions for processing observational data
- "visualizations.ipynb": Jupyter notebook, which produces all the figures of the paper (see python_env.yaml for the conda environment used)