#' ---
#' title: ABCDCon Control Analysis for Overly-influential Voxels
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
project_dir <- "../"
config <- yaml::yaml.load_file(paste0(project_dir,"config.yml"))

#' ## Set paths as variables
analyzed_mri_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_mri))
raw_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$raw_behavioral))
analyzed_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_behavioral))
roi_dirs <- c("ashs_left","ashs_right")
num_vox <- 5 # Set the number of voxels to pull

#' ## Setup other variables
subjects <- c(1:30)
exclude_subjects <- c(3:5,15:17)
subjects <- subjects[!is.element(subjects, exclude_subjects)]

rois <- c("brCA1_body.nii", "brCA2_3_DG_body.nii")

#' ### Flags
GRAPH_FLAG <- 0

#' # Loop and load in the data
#+ warning = FALSE
for(idir in c(1:length(roi_dirs))) {

  for(iroi in c(1:length(rois))){

    for(isub in subjects) {

      # define subject specific directories and filenames
      subj_number <- halle::format_three_digit_subject_id(isub)
      base_roi_dir <- paste(analyzed_mri_dir,halle::ensure_trailing_slash(subj_number),halle::ensure_trailing_slash('ROIs'),sep="")
      cur_roi_dir <- paste(base_roi_dir,halle::ensure_trailing_slash(roi_dirs[idir]),sep="")
      cur_roi <- strsplit(rois[iroi],".nii")

      # figure out L vs R if using ASHS ROIs
      if(roi_dirs[idir]=="ashs_left"){
        hemi_label <- "left"
      } else if(roi_dirs[idir]=="ashs_right") {
        hemi_label <- "right"
      }

      # test to make sure the current file exists
      # if not, skip and continue on
      # this file is created by `RSA_btwn_runs_exclude_outlier_trials.m`
      cur_file <- paste(
        cur_roi_dir, cur_roi,'_pattern_mtx_no_outlier_trials_all_runs.mat',sep="")

      if(file.exists(cur_file)){

        # --- read in the subject's data for the current ROI ---
        subj_pattern_mtx <- data.matrix(as.data.frame(R.matlab::readMat(cur_file)))

        # --- take abs value ---
        #  we assume that direction is essentially meaningless when not compared to a baseline
        abs_mtx <- subj_pattern_mtx %>%
          abs()

        # --- assume stationary across time ---
        #  if we are concerned about just a few voxels driving a pattern, this assumption seems reasonable
        mean_mtx <- abs_mtx %>%
          as.data.frame() %>%
          dplyr::mutate(mean_vox_val = rowMeans(., na.rm = TRUE)) %>%
          # add a column with voxel number
          # based on: https://stackoverflow.com/questions/23518605/add-an-index-numeric-id-column-to-large-data-frame-in-r
          tibble::rowid_to_column(., "voxel_number")

        # --- rank order voxels by mean value ---
        # test this ( remember that default is using the last column):
        # df <- data.frame(x = c(10, 4, 1, 6, 3, 1, 1), y = c(1, 2, 4, 6, 8, 4, 5), z = c(10, 9, 8, 4, 6, 5, 4)) %>% tibble::rowid_to_column("id")
        # df %>%
        #  dplyr::top_n(4)

        rank_vals <- mean_mtx %>%
          dplyr::top_n(n = num_vox, wt = mean_vox_val)

        exclude_vox_ids <- rank_vals$voxel_number

        # ---   save out voxels to exclude ---
        write.table(exclude_vox_ids, file = file.path(analyzed_mri_dir, subj_number, sprintf("%s_%s_top_%s_voxels.csv", hemi_label, cur_roi, num_vox)), row.names=FALSE, col.names = FALSE, sep = ",")

      } else {
        print(sprintf("Current pattern matrix file %s does not exist. Continuing on.", cur_file))
        next
      }#if cur_file exists

    } #isub

  } #iroi

} #idir
