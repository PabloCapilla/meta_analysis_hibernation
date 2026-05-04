# Data and code repository for:

**Effects of animal dormancy on oxidative stress, immune status, and glucocorticoids: a meta-analysis**

Pablo Burraco, Tamara G Petrović, Pablo Capilla-Lasheras, Marko Prokić. [Link to bioRxiv](https://www.biorxiv.org/content/10.1101/2025.10.02.680078v1.full.pdf)

# Repositoty structure:

- **01_data**:

  - *processed_RDS_data_files*: RDS files produced by [prep'ing script](LINK). These are the data files used in subsequent analyses.

  - *raw_data*: raw datatables to start analysis.

- **02_scripts**:

  - *01_data_prep*: scripts to process (raw data)\[LINK\].

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

- **05_reports**: html files showing content of each (analysis script)\[LINK\].