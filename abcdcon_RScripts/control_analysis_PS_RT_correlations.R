#' ---
#' title: ABCDCon RT Control Analysis
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
graph_fpath_out <- paste0(halle::ensure_trailing_slash(dropbox_dir),
                          halle::ensure_trailing_slash("writeups"),
                          halle::ensure_trailing_slash("figures"))

#' ## Setup other variables
#' ### Flags
SAVE_GRAPHS_FLAG <-1

#' # Load in PS data
load(file.path(analyzed_mri_dir, 'group_z_renamed_spatial_temporal_PS_by_trial.RData'))

#' ## Peek at the PS data
head(all_trials_z_better_names)
unique(all_trials_z_better_names$roi)
unique(all_trials_z_better_names$condition)
unique(all_trials_z_better_names$subj)

#' # Load behavioral data
load(paste0(raw_behavioral_dir,'group_data.RData'))

#' ## Rename it so not called `data`
behav_dat <- data %>%
  # format subject IDs so can merge w/ PS data
  dplyr::mutate(subj = sprintf('s%03d', subj_num))

#' ## Take a peek
head(behav_dat)
unique(behav_dat$subj)

#' # Merge behavioral data and PS data using `subj`
alldat <- dplyr::full_join(all_trials_z_better_names, behav_dat, by = "subj")

#' ##
