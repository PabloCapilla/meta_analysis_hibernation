![https://zenodo.org/badge/doi/10.5281/zenodo.22255960.svg](https://zenodo.org/badge/doi/10.5281/zenodo.22255960.svg)




# Data and code repository for:

**Effects of animal dormancy on oxidative stress, immune status, and glucocorticoids: a meta-analysis**

Pablo Burraco, Tamara G Petrović, Pablo Capilla-Lasheras, Marko Prokić. [![DOI:10.1101/2025.10.02.680078v1.full](http://img.shields.io/badge/DOI-10.1101/2025.10.02.680078v1.full.svg)](https://doi.org/10.1101/2025.10.02.680078v1.full) 



# Repositoty structure:

- **01_data**:

  - *processed_RDS_data_files*: RDS files produced by [prep'ing script](https://github.com/PabloCapilla/meta_analysis_hibernation/tree/main/02_scripts/01_data_prep). These are the data files used in subsequent analyses.

  - *raw_data*: raw datatables to start analysis.

- **02_scripts**:

  - *01_data_prep*: scripts to process [raw data](https://github.com/PabloCapilla/meta_analysis_hibernation/tree/main/01_data/raw_data).

  - *02_phylo_cor*: scripts to generate phylogenetic var-covar matrices for phylogenetically-controlled meta-analysis.

  - *03_analysis*: R markdown documents with all models and plots.

- **03_plots**:

  - *01_hibernation*

  - *02_arousal*

  - *03_aestivation*

- **04_models**:

  - *01_hibernation*

  - *02_arousal*

  - *03_aestivation*

- **05_reports**: html files showing content of each [analysis script](https://github.com/PabloCapilla/meta_analysis_hibernation/tree/main/02_scripts/03_analysis)