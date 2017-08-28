#' ---
#' title: ABCDCon control analysis plot PS correlations
#' author: Halle R. Dimsdale-Zucker
#' output:
#'  html_document:
#'    toc: true
#'    toc_depth: 4
#'    toc_float:
#'      collapsed: false
#'      smooth_scroll: false
#'    number_sections: true
#'    theme: spacelab
#' ---

#+ initialize, warning = FALSE, message = FALSE
devtools::load_all()

# add dplyr with library()
# NB: this is non-standard
# correct way would be to make this script into a function,
# store in the R/ directory,
# and use @importFrom dplyr "%>%",
# but this is incompatible with getting in-line results with code in knitr output
library(dplyr)

#' # Setup
#' ## Load config file
project_dir <- ("../")
config <- yaml::yaml.load_file(paste0(project_dir, "config.yml"))

#' ## Set paths as variables
analyzed_mri_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_mri))
raw_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$raw_behavioral))
analyzed_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_behavioral))
dropbox_dir <- halle::ensure_trailing_slash(config$directories$dropbox_abcdcon)
dropbox_graph_fpath_out <- paste0(halle::ensure_trailing_slash(dropbox_dir),
                                  halle::ensure_trailing_slash("writeups"),
                                  halle::ensure_trailing_slash("figures"))
graph_fpath_out <- paste0(dropbox_graph_fpath_out, "plot-PS-matrices")
dir.create(graph_fpath_out) # will throw an error if this already exists

#' ## Setup other variables
#' ### Flags
SAVE_GRAPHS_FLAG <-1

#+ label="Load data"
#' ## Load in PS data
# this file is saved out in `mixed_models.R`
load(paste0(analyzed_mri_dir, 'group_z_renamed_spatial_temporal_PS_by_trial.RData'))

#' # Filter out data that's not needed
tidy_trials <- all_trials_z_better_names %>%
  dplyr::filter(hemi == "left",
                condition %in% c("diffVideo_diffHouse", "sameVideo_sameHouse", "diffVideo_sameHouse"),
                roi %in% c("CA1_body", "CA2_3_DG_body")) %>%
  dplyr::select(-z_r)

#' # Loop and create plots
subjects <- unique(tidy_trials$subj)
nsub <- length(subjects)
conditions <- unique(tidy_trials$condition)
ncond <- length(conditions)
rois <- unique(tidy_trials$roi)
nroi <- length(rois)

for(isubj in 1:nsub){
  for(icond in 1:ncond){
    for(iroi in 1:nroi){

      cur_subj <- subjects[isubj]
      cur_cond <- conditions[icond]
      cur_roi <- rois[iroi]

      cur_dat <- tidy_trials %>%
        dplyr::filter(subj == cur_subj,
                      condition == cur_cond,
                      roi == cur_roi)

      cur_dat_fmt <- cur_dat %>%
        tidyr::spread(col_name, r) %>%
        dplyr::select(-subj, -roi, -hemi, -condition, -row_name)

      # TODO: format plots so subjects are in columns and conditions are in rows. save out separate files for each ROI
      GGally::ggcorr(cur_dat_fmt, size = 0)
    } #iroi
  } #icond
} #isubj
