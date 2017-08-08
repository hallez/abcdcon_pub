#' ---
#' title: ABCDCon Beta Distribution Control Analysis
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
SAVE_GRAPHS_FLAG <- 1

#' ### Included subjects
subjects <- c(1:2, 6:14, 18:21, 23:30)
subj_formatted <- NULL
for(isub in 1:length(subjects)){
  subj_formatted[isub] <- halle::format_three_digit_subject_id(subjects[isub])
}

#' ### Other variables (e.g., ROIs of interest)
hemis <- c("ashs_left", "ashs_right")
ROIs <- c("brCA1_body", "brCA2_3_DG_body")

#' # Load in data
cur_betas_fpath <- file.path(graph_fpath_out, 'allSubj_allROI_means.mat')
if(file.exists(cur_betas_fpath)){
  cur_betas <- as.data.frame(R.matlab::readMat(cur_betas_fpath))

  # --- Label `cur_betas` w/ subject IDs ---
  # based on: https://stackoverflow.com/questions/22670541/subsetting-a-matrix-by-row-names
  cur_ids <- cur_betas[row.names(cur_betas) == "ids", ]
  cur_ids_fmt <- cur_ids %>%
    tidyr::gather(colname_junk, subj_id) %>%
    dplyr::mutate(subj_id_unlist = unlist(subj_id))

  cur_betas_w_ids <- cur_betas
  colnames(cur_betas_w_ids) <- cur_ids_fmt$subj_id_unlist

  # --- Flip dataframe around so have rows instead of cols ---
  cur_betas_tidy <- cur_betas_w_ids %>%
    dplyr::mutate(lbls = rownames(.)) %>%
    tidyr::gather(subj_id, val, -lbls) %>%
    # FIGURE OUT HOW TO BREAK DOWN `LBLS` INTO COLUMN NAMES (e.g., eventually want to have columns like `hemi`, `roi`, `value_type`)
    tidyr::spread(key = lbls, value = val) %>%
    dplyr::select(-ids)




  } else {
    print(sprintf("Current mean betas file %s does not exist. Continuing on.", cur_betas_fpath))
    next
} #if(file.exists


#' # Tidy up group dataframe

#' # Plot

