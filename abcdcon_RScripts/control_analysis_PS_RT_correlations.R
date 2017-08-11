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

#' ## Split out trials from trial pairs so can merge w/ behavioral data
PS_trial_ids <- all_trials_z_better_names %>%
  # NB: will get a warning that there are `too many values`, but still separates correctly (I think this is b/c there are multiple fiels for `sep`)
  tidyr::separate(row_name, into = c("row.run_trialnum", "row.trial_lbl"), sep = "Brown|Gray", remove = FALSE) %>%
  dplyr::mutate(row.runID_trialID = gsub("(^Run)_([01234]*)_(Trial)_([0123456789]*)_", "\\1\\2_\\3\\4", row.run_trialnum)) %>%
  tidyr::separate(col_name, into = c("col.run_trialnum", "col.trial_lbl"), sep = "Brown|Gray", remove = FALSE) %>%
  dplyr::mutate(col.runID_trialID = gsub("(^Run)_([01234]*)_(Trial)_([0123456789]*)_", "\\1\\2_\\3\\4", col.run_trialnum)) %>%
  dplyr::select(-row.trial_lbl, -row.run_trialnum, -col.trial_lbl, -col.run_trialnum)

#' # Load behavioral data
load(paste0(raw_behavioral_dir,'group_data.RData'))

#' ## Rename it so not called `data`
behav_dat <- data %>%
  # format subject IDs so can merge w/ PS data
  dplyr::mutate(subj = sprintf('s%03d', subj_num))

#' ## Check the data before messing w/ it
head(behav_dat)

#' ## Format so can merge w/ PS data
behav_trim <- behav_dat %>%
  # create a column that combines run and trial number so can merge w/ PS data
  # duplicate this column so can match w/ either row or col is PS data
  dplyr::mutate(row.runID_trialID = sprintf('Run%02s_Trial%03s', BlockID, objRec_trialNum)) %>%
  dplyr::mutate(col.runID_trialID = sprintf('Run%02s_Trial%03s', BlockID, objRec_trialNum)) %>%
  # only select RT column
  dplyr::select(subj, col.runID_trialID, row.runID_trialID, MemRT)

#' ## Take a peek
head(behav_trim)

#' # Merge behavioral data and PS data
alldat <- dplyr::left_join(PS_trial_ids, behav_trim, by = c("subj", "row.runID_trialID")) %>%
  # rename colID column b/c it gets messed up in join
  dplyr::rename(col.runID_trialID = col.runID_trialID.x) %>%
  # rename MemRT column so it doesn't get crushed when do left join again
  dplyr::rename(MemRT_row = MemRT) %>%
  dplyr::left_join(., behav_trim, by = c("subj", "col.runID_trialID")) %>%
  dplyr::rename(MemRT_col = MemRT,
                row.runID_trialID = row.runID_trialID.x) %>%
  # remove columns that get duplicated
  dplyr::select(-col.runID_trialID.y, -row.runID_trialID.y) %>%
  # compute RT difference for each trial pair
  dplyr::mutate(RT_diff = MemRT_row - MemRT_col)

#' ## Sanity spot check RTs
#' ### Rows
alldat %>%
  dplyr::filter(row.runID_trialID == "Run02_Trial010",
                hemi == "left",
                roi == "CA1",
                subj == "s006") %>%
  dplyr::distinct(MemRT_row)

behav_trim %>%
  dplyr::filter(row.runID_trialID == "Run02_Trial010",
                subj == "s006") %>%
  dplyr::distinct(MemRT)

#' ### Columns
alldat %>%
  dplyr::filter(col.runID_trialID == "Run04_Trial038",
                hemi == "left",
                roi == "CA1",
                subj == "s006") %>%
  dplyr::distinct(MemRT_col)

behav_trim %>%
  dplyr::filter(col.runID_trialID == "Run04_Trial038",
                subj == "s006") %>%
  dplyr::distinct(MemRT)

#' # Plot PS values against RT difference
#' ## Define functions
scatter_cond <- function(dat_in){
  p <- NULL
  p <<- dat_in %>%
    # only take positive RT differences b/c the data are symmetric (b/c trial pairs are all rows and all columns)
    dplyr::filter(RT_diff >=0) %>%
    ggplot2::ggplot(ggplot2::aes(x = r, y = RT_diff)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = lm)
}

scatter_cond_facet_by_subj <- function(dat_in){
  p <- NULL
  p <<- dat_in %>%
    # only take positive RT differences b/c the data are symmetric (b/c trial pairs are all rows and all columns)
    dplyr::filter(RT_diff >=0) %>%
    ggplot2::ggplot(ggplot2::aes(x = r, y = RT_diff)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~subj) +
    ggplot2::geom_smooth(method = lm)
}

#' ## Plot each condition
# TODO: Create a for loop and autogenerate filenames
alldat %>%
  dplyr::filter(hemi == "left",
                roi == "CA1",
                condition == "sameVideo_sameHouse") %>%
  scatter_cond()

print(p)

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = paste0(graph_fpath_out,
                                "RT-PS_correlation_sameVideo_sameHouse.pdf"),
                  width=8, height=6)
}

alldat %>%
  dplyr::filter(hemi == "left",
                roi == "CA1",
                condition == "sameVideo_sameHouse") %>%
  scatter_cond_facet_by_subj()

print(p)

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = paste0(graph_fpath_out,
                                "RT-PS_correlation_sameVideo_sameHouse_bySubj.pdf"),
                  width=8, height=6)
}



