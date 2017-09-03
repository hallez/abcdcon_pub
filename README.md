ABCDCon
===

# Repository details
This repository contains the code used to collect data and generate analyses for the following paper:

[Dimsdale-Zucker, H. R., Ritchey, M., Ekstrom, A. D., Yonelinas, A. P., Ranganath, C. (2017). CA1 and CA3 differentially support spontaneous retrieval of episodic contexts within human hippocampal subfields. bioRxiv. doi: https://doi.org/10.1101/142349](http://www.biorxiv.org/content/early/2017/05/25/142349)

The scripts contained here were used to collect and analyze data using Matlab r2014b and R version 3.3.2 run in RStudio and have not been tested for compatibility with later versions or across operating systems. They are provided for the purposes of openness and transparency of the data generation and analysis process.

**Note, these are the scripts I used *at the time* of data collection and analysis and I might suggest more efficient or cleaner coding solutions to the same problems now (so use at your own risk and please be kind if you find bugs).**

# Environment setup
## Setting paths
1. Copy the config.yml.example file to config.yml and edit the paths for the current machine.
```
cp config.yml.example config.yml
vim config.yml
```

## Python
1. Install python
2. Install python requirements:
```
pip install -r requirements.txt
```

## R
* Required pacakges:
  * dplyr (>= 0.4.0)
  * ggplot2 (>= 1.0.1)
  * [halle (>=0.5.6)](https://github.com/hallez/halle)
  * lme4(>= 1.1-11)
  * R.matlab (>= 3.2.0)
  * reshape2(>= 1.4.1)
  * stringr (>= 1.0.0)
  * tidyr (>= 0.2.0)
  * yaml (>= 2.1.13)

## Outside functions and packages
Code from other sources has been included in the `vendor` directory. This is not code that I wrote and therefore is subject to the usage license and instructions of those authors. It is included here out of convenience for those trying to run the scripts in the current repository.

# Scripts
Variables that are common across Matlab scripts will be set by `scripts/mri_analyses/initialize_ABCDCon.m` (This will need to be added to your Matlab path).
## Data collection
* These scripts are all contained in `scripts/run_task`
* For the study, scripts were run in the following order:
  1. `ABCDCon_contextEnc.m`
  2. `ABCDCon_objectEnc.m`
  3. `ABCDCon_objectRecog_MRI.m`
  4. `ABCDCon_locationRecog.m`

## Behavioral data analysis
These scripts can be found in `abcdcon_RScripts`. They are setup as an [R package](http://r-pkgs.had.co.nz/) for ease of loading required packages and functions, but this is probably a non-standard use of an R package. You should open `abcdcon.Rproj` and then run scripts from here (this should set the current working directory correctly without using `setwd()`).
1. `load_ABCDCon.R`
2. `analyze_ABCDCon.R`

## MRI: Univariate analyses
Assumes have already run `load_ABCDCon.R`, downloaded .zip MRI files, converted dicom images, set subject-specific folder names, commonized folder names (e.g., `run1` instead of whatever the MRI scanner outputs), and run [data quality assurance](https://github.com/ritcheym/DML_QA). Scripts assume they are being run relative to `scripts/mri_analyses`.
1. Preprocess the data using `preprocess.m`
2. Generate regressors for the first-level analyses using `generate_regressors.m`
3. Run first-level analyses using `first_level.m`
4. Put contrast images into MNI space using `new_segment_write_deformations_batch.m`
5. Smooth the resultant `wcon` images using a 3mm kernal. (I set this up manually in the SPM GUI.)
6. Run second-level (group) analyses using `second_level_job.m`

## MRI: Multivariate analyses 
Assumes you have already preprocessed, run QA (to generate spike regressors), and have traced ROIs.
1. Generate single-trial regressors using `RSA_generate_single_trial_regressors.m`
2. Estimate the single trial betas using `RSA_single_trial_models_batch.m`
3. Identify outlier beta timepoints using `RSA_beta_timeseries_graphs.m` to visually identify the group threshold and then `RSA_beta_timeseries_id_outliers.m` to mark excluded betas on a subject-by-subject basis
3. Gather the ROIs of interest ensuring that have them split so can look at body separate from head and tail. This involves extracting the ROIs of interest from the ASHS tracing files (`RSA_extractROIs.m`), combining subfields (`RSA_combine_rois.py`), and splitting along the long axis (`RSA_split_anterior_posterior.m`).
   * In order to split into head/body/tail, you must manually identify transitions. The slice numbers that define the subject-unique transitions are written out into files; the splitting scripts will assume the transitions file is in the <RAW_BEHAVIORAL> directory for each subject and are called `s*_hc_transitions.yml`. Use the `hc_transitions.yml.example` file for consistent formatting.
4. Reslice the ROIs into EPI space (if using ASHS, these will be in T2 space) using `RSA_reslice_t2_and_ROIs_batch.m`
5. Binarize the ROIs so they can be used as masks using `RSA_binarize_ROIs_batch.m`
6. Generate trial labels
   * This is a convoluted system that should be re-implemented differently if these analysis scripts were re-generated from scratch. However, what I did was first extract trial IDs from the single trial betas using `RSA_trial_ids_from_betas.m` and then use `decode_single_trial_labels.R` to mark trials that needed to be excluded (e.g., outlier betas).
6. Extract beta values for trials of interest using `RSA_btwn_runs_exclude_outlier_trials.m`
8. Calculate the pattern similarity values of interest using `pattern_similarity_no_outlier_trials_load_data_btwn_runs.R`
9. Analyze using mixed models with `mixed_models.R`

## Supplemental analyses 
1. Use FIR to estimate shape of univariate response in each subfield (Supplemental Figure X):
  1. (matlab scripts here)
  1. Plot the data: `control_analysis_FIR_to_GLM_plot_betas.R`
1. Look at distribution of beta values (Supplemental Figure X)
  1. Get trial-wise mean betas: `control_analysis_betas_by_ROI.m`
  1. Plot distrubtions: `control_analysis_plot_betas_by_ROI.m`
1. Randomly select trials so that all conditions have equal bin sizes: `control_analysis_matched_trial_numbers.R`
1. Remove voxels of influence:
  1. Identify influential voxels: `control_analysis_drop_voxels.R`
  1. Remove these voxels from pattern matrices: `control_analysis_remove_top_voxels.m`
  1. Re-compute PS without these voxels: `control_analysis_PS_truncated_voxels.R`
  1. Run stats: `control_analysis_PS_truncated_voxels_mixed_models.R`

