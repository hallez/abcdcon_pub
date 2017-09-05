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

        # OPTION 1: select largest voxels
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
        write.table(exclude_vox_ids, file = file.path(analyzed_mri_dir, subj_number, sprintf("%s_%s_top_%s_voxels_by_mean.csv", hemi_label, cur_roi, num_vox)), row.names=FALSE, col.names = FALSE, sep = ",")

        # OPTION 2: select most variable voxels
        # --- find SD of each voxel across time ---
        sd_mtx <- subj_pattern_mtx %>%
          as.data.frame() %>%
          tibble::rowid_to_column(., "voxel_number") %>%
          # have multiple rows per voxel (ie, long data instead of wide)
          tidyr::gather(key = timepoint_id, value = timepoint_value, -voxel_number) %>%
          # group by voxel ID
          dplyr::group_by(voxel_number) %>%
          # take SD for each voxel
          dplyr::summarise(sd_vox_val = sd(timepoint_value, na.rm = TRUE))

        # --- rank by SD ---
        rank_sds <- sd_mtx %>%
          dplyr::top_n(n = num_vox, wt = sd_vox_val)

        exclude_sd_vox_ids <- rank_sds$voxel_number

        # --- save out voxels to exclude ---
        write.table(exclude_sd_vox_ids, file = file.path(analyzed_mri_dir, subj_number, sprintf("%s_%s_top_%s_voxels_by_SD.csv", hemi_label, cur_roi, num_vox)), row.names=FALSE, col.names = FALSE, sep = ",")

        # --- OPTION #3 ---
        # take the voxels that are the most variable between conditions (same video same house vs. different video same house)
        # based on: https://stackoverflow.com/questions/22446825/perform-pairwise-comparison-of-matrix
        # try it out w/ a subset of the data to ensure this procedure works:
        # d <- subj_pattern_mtx[1:3,1:5]
        # d_diff <- apply(combn(ncol(d), 2), 2, function(x) (d[,x[1]] - d[,x[2]])^2)
        # colnames(d_diff) <- apply(combn(ncol(d), 2), 2, function(x) paste(colnames(d)[x], collapse=' - '))

        # index based on trial pair ID
        # these tidied ID files are created by `control_analysis_plot_PS_matrices.R`
        tidy_ids_file <- file.path(analyzed_mri_dir, subj_number, "ROIs",
                                   "ashs_left", sprintf('%s_pattern_mtx_ids_tidied_no_outlier_trials_all_runs.RData',cur_roi))
        if(file.exists(tidy_ids_file)){
          load(tidy_ids_file)

          # --- relabel patterns by video ID ---
          # this would be easier had I used matrices to index patterns, but I didn't so this is hacky
          subj_pattern_mtx_by_video <- subj_pattern_mtx
          colnames(subj_pattern_mtx_by_video) <- subj_pattern_ids_tidy$video_id

          # --- reset NaN values to 0 ---
          # is this the appropriate way to handle NaNs in this case?
          # if don't reset to 0, then all values in resultant matrix are NAs
          subj_pattern_mtx_by_video_no_NaN <- subj_pattern_mtx_by_video
          subj_pattern_mtx_by_video_no_NaN[is.na(subj_pattern_mtx_by_video_no_NaN)] <- 0

          # --- compute squared differences between every trial ---
          # could just take differences, but squaring makes sign invariant
          pattern_diffs_by_video_no_NaN <- apply(combn(ncol(subj_pattern_mtx_by_video_no_NaN), 2), 2, function(x) (subj_pattern_mtx_by_video_no_NaN[,x[1]] - subj_pattern_mtx_by_video_no_NaN[,x[2]])^2)
          # also put on column names that reflect which trials are being subtracted and squared
          # this is critical b/c it enables the hacky filtering of trials
          colnames(pattern_diffs_by_video_no_NaN) <- apply(combn(ncol(subj_pattern_mtx_by_video_no_NaN), 2), 2, function(x) paste(colnames(subj_pattern_mtx_by_video_no_NaN)[x], collapse='.'))

          # spot check a few
          subj_pattern_mtx_by_video_no_NaN[1:3,1:5]
          pattern_diffs_by_video_no_NaN[1:3,1:5]

          # --- grab same video, same house trials ---
          SVSH_names <- c("Brown1.Brown1", "Brown2.Brown2", "Brown3.Brown3", "Brown4.Brown4", "Brown5.Brown5", "Brown6.Brown6",
                          "Brown7.Brown7", "Brown8.Brown8", "Brown9.Brown9", "Brown10.Brown10", "Brown11.Brown11", "Brown12.Brown12",
                          "Gray1.Gray1", "Gray2.Gray2", "Gray3.Gray3", "Gray4.Gray4", "Gray5.Gray5", "Gray6.Gray6",
                          "Gray7.Gray7", "Gray8.Gray8", "Gray9.Gray9", "Gray10.Gray10", "Gray11.Gray11", "Gray12.Gray12")
          # get column IDs that match these possible names
          SVSH_col_ids <- which(is.element(colnames(pattern_diffs_by_video_no_NaN), SVSH_names))
          # now filter so just have these columns
          pattern_diffs_by_video_SVSH <- pattern_diffs_by_video_no_NaN[,SVSH_col_ids]

          # --- grab indices of different video, same house trials ---
          # even though brown1-brown2 should be the same as brown2-brown1, grab both
          DVSH_names <- c("Brown1.Brown2", "Brown1.Brown3", "Brown1.Brown4", "Brown1.Brown5", "Brown1.Brown6",
                          "Brown1.Brown7", "Brown1.Brown8", "Brown1.Brown9", "Brown1.Brown10", "Brown1.Brown11", "Brown1.Brown12",
                          "Brown2.Brown1", "Brown2.Brown3", "Brown2.Brown4", "Brown2.Brown5", "Brown2.Brown6",
                          "Brown2.Brown7", "Brown2.Brown8", "Brown2.Brown9", "Brown2.Brown10", "Brown2.Brown11", "Brown2.Brown12",
                          "Brown3.Brown1", "Brown3.Brown2", "Brown3.Brown4", "Brown3.Brown5", "Brown3.Brown6",
                          "Brown3.Brown7", "Brown3.Brown8", "Brown3.Brown9", "Brown3.Brown10", "Brown3.Brown11", "Brown3.Brown12",
                          "Brown4.Brown1", "Brown4.Brown2", "Brown4.Brown3", "Brown4.Brown5", "Brown4.Brown6",
                          "Brown4.Brown7", "Brown4.Brown8", "Brown4.Brown9", "Brown4.Brown10", "Brown4.Brown11", "Brown4.Brown12",
                          "Brown5.Brown1", "Brown5.Brown2", "Brown5.Brown3", "Brown5.Brown4", "Brown5.Brown6",
                          "Brown5.Brown7", "Brown5.Brown8", "Brown5.Brown9", "Brown5.Brown10", "Brown5.Brown11", "Brown5.Brown12",
                          "Brown6.Brown1", "Brown6.Brown2", "Brown6.Brown3", "Brown6.Brown4", "Brown6.Brown5",
                          "Brown6.Brown7", "Brown6.Brown8", "Brown6.Brown9", "Brown6.Brown10", "Brown6.Brown11", "Brown6.Brown12",
                          "Brown7.Brown1", "Brown7.Brown2", "Brown7.Brown3", "Brown7.Brown4", "Brown7.Brown5", "Brown7.Brown6",
                          "Brown7.Brown8", "Brown7.Brown9", "Brown7.Brown10", "Brown7.Brown11", "Brown7.Brown12",
                          "Brown8.Brown1", "Brown8.Brown2", "Brown8.Brown3", "Brown8.Brown4", "Brown8.Brown5", "Brown8.Brown6",
                          "Brown8.Brown7", "Brown8.Brown9", "Brown8.Brown10", "Brown8.Brown11", "Brown8.Brown12",
                          "Brown9.Brown1", "Brown9.Brown2", "Brown9.Brown3", "Brown9.Brown4", "Brown9.Brown5", "Brown9.Brown6",
                          "Brown9.Brown7", "Brown9.Brown8", "Brown9.Brown10", "Brown9.Brown11", "Brown9.Brown12",
                          "Brown10.Brown1", "Brown10.Brown2", "Brown10.Brown3", "Brown10.Brown4", "Brown10.Brown5", "Brown10.Brown6",
                          "Brown10.Brown7", "Brown10.Brown8", "Brown10.Brown9", "Brown10.Brown11", "Brown10.Brown12",
                          "Brown11.Brown1", "Brown11.Brown2", "Brown11.Brown3", "Brown11.Brown4", "Brown11.Brown5", "Brown11.Brown6",
                          "Brown11.Brown7", "Brown11.Brown8", "Brown11.Brown9", "Brown11.Brown10", "Brown11.Brown12",
                          "Brown12.Brown1", "Brown12.Brown2", "Brown12.Brown3", "Brown12.Brown4", "Brown12.Brown5", "Brown12.Brown6",
                          "Brown12.Brown7", "Brown12.Brown8", "Brown12.Brown9", "Brown12.Brown10", "Brown12.Brown11",
                          "Gray1.Gray2", "Gray1.Gray3", "Gray1.Gray4", "Gray1.Gray5", "Gray1.Gray6",
                          "Gray1.Gray7", "Gray1.Gray8", "Gray1.Gray9", "Gray1.Gray10", "Gray1.Gray11", "Gray1.Gray12",
                          "Gray2.Gray1", "Gray2.Gray3", "Gray2.Gray4", "Gray2.Gray5", "Gray2.Gray6",
                          "Gray2.Gray7", "Gray2.Gray8", "Gray2.Gray9", "Gray2.Gray10", "Gray2.Gray11", "Gray2.Gray12",
                          "Gray3.Gray1", "Gray3.Gray2", "Gray3.Gray4", "Gray3.Gray5", "Gray3.Gray6",
                          "Gray3.Gray7", "Gray3.Gray8", "Gray3.Gray9", "Gray3.Gray10", "Gray3.Gray11", "Gray3.Gray12",
                          "Gray4.Gray1", "Gray4.Gray2", "Gray4.Gray3", "Gray4.Gray5", "Gray4.Gray6",
                          "Gray4.Gray7", "Gray4.Gray8", "Gray4.Gray9", "Gray4.Gray10", "Gray4.Gray11", "Gray4.Gray12",
                          "Gray5.Gray1", "Gray5.Gray2", "Gray5.Gray3", "Gray5.Gray4", "Gray5.Gray6",
                          "Gray5.Gray7", "Gray5.Gray8", "Gray5.Gray9", "Gray5.Gray10", "Gray5.Gray11", "Gray5.Gray12",
                          "Gray6.Gray1", "Gray6.Gray2", "Gray6.Gray3", "Gray6.Gray4", "Gray6.Gray5",
                          "Gray6.Gray7", "Gray6.Gray8", "Gray6.Gray9", "Gray6.Gray10", "Gray6.Gray11", "Gray6.Gray12",
                          "Gray7.Gray1", "Gray7.Gray2", "Gray7.Gray3", "Gray7.Gray4", "Gray7.Gray5", "Gray7.Gray6",
                          "Gray7.Gray8", "Gray7.Gray9", "Gray7.Gray10", "Gray7.Gray11", "Gray7.Gray12",
                          "Gray8.Gray1", "Gray8.Gray2", "Gray8.Gray3", "Gray8.Gray4", "Gray8.Gray5", "Gray8.Gray6",
                          "Gray8.Gray7", "Gray8.Gray9", "Gray8.Gray10", "Gray8.Gray11", "Gray8.Gray12",
                          "Gray9.Gray1", "Gray9.Gray2", "Gray9.Gray3", "Gray9.Gray4", "Gray9.Gray5", "Gray9.Gray6",
                          "Gray9.Gray7", "Gray9.Gray8", "Gray9.Gray10", "Gray9.Gray11", "Gray9.Gray12",
                          "Gray10.Gray1", "Gray10.Gray2", "Gray10.Gray3", "Gray10.Gray4", "Gray10.Gray5", "Gray10.Gray6",
                          "Gray10.Gray7", "Gray10.Gray8", "Gray10.Gray9", "Gray10.Gray11", "Gray10.Gray12",
                          "Gray11.Gray1", "Gray11.Gray2", "Gray11.Gray3", "Gray11.Gray4", "Gray11.Gray5", "Gray11.Gray6",
                          "Gray11.Gray7", "Gray11.Gray8", "Gray11.Gray9", "Gray11.Gray10", "Gray11.Gray12",
                          "Gray12.Gray1", "Gray12.Gray2", "Gray12.Gray3", "Gray12.Gray4", "Gray12.Gray5", "Gray12.Gray6",
                          "Gray12.Gray7", "Gray12.Gray8", "Gray12.Gray9", "Gray12.Gray10", "Gray12.Gray11")
          DVSH_col_ids <- which(is.element(colnames(pattern_diffs_by_video_no_NaN), DVSH_names))
          pattern_diffs_by_video_DVSH <- pattern_diffs_by_video_no_NaN[,DVSH_col_ids]

          # --- figure out the voxels with the largest squared difference between conditions ---
          # need to remove column names b/c it doesn't like that they're duplicated
          pattern_diffs_by_video_SVSH_no_col_ids <- pattern_diffs_by_video_SVSH
          colnames(pattern_diffs_by_video_SVSH_no_col_ids) <- NULL

          SVSH_means <- pattern_diffs_by_video_SVSH_no_col_ids %>%
            as.data.frame() %>%
            # compute a mean value for each voxel w/in a condition
            dplyr::mutate(mean_vox_val_SVSH = rowMeans(., na.rm = TRUE)) %>%
            tibble::rowid_to_column(., "voxel_number") %>%
            dplyr::select(voxel_number, mean_vox_val_SVSH)

          # repeat those steps for the DVSH condition
          pattern_diffs_by_video_DVSH_no_col_ids <- pattern_diffs_by_video_DVSH
          colnames(pattern_diffs_by_video_DVSH_no_col_ids) <- NULL

          DVSH_means <- pattern_diffs_by_video_DVSH_no_col_ids %>%
            as.data.frame() %>%
            dplyr::mutate(mean_vox_val_DVSH = rowMeans(., na.rm = TRUE)) %>%
            tibble::rowid_to_column(., "voxel_number") %>%
            dplyr::select(voxel_number, mean_vox_val_DVSH)

          # merge the dataframes
          SVSH_DVSH <- SVSH_means %>%
            dplyr::left_join(DVSH_means) %>%
            # compute the difference in each voxel's mean squared difference value between conditions
            dplyr::mutate(vox_diff = mean_vox_val_SVSH - mean_vox_val_DVSH,
                          vox_sq_diff = vox_diff^2)

          # rank voxels by largest mean squared difference between conditions
          rank_sq_condition_diffs <- SVSH_DVSH %>%
            dplyr::top_n(n = num_vox, wt = vox_sq_diff)

          exclude_cond_diffs_vox_ids <- rank_sq_condition_diffs$voxel_number

          # --- save out voxels to exclude ---
          write.table(exclude_cond_diffs_vox_ids, file = file.path(analyzed_mri_dir, subj_number, sprintf("%s_%s_top_%s_voxels_by_condition_diff.csv", hemi_label, cur_roi, num_vox)), row.names=FALSE, col.names = FALSE, sep = ",")

        } else {
          print(sprintf("Current tidied ids file %s does not exist. Continuing on.", tidy_ids_file))
          next
        }#file.exists(tidy_ids_file)

        # TODO: FIGURE OUT WAY TO QUANTIFY OVERLAP IN THESE APPROACHES - would be interesting to see which voxels are the same ones that get removed in each

      } else {
        print(sprintf("Current pattern matrix file %s does not exist. Continuing on.", cur_file))
        next
      }#if cur_file exists

    } #isub

  } #iroi

} #idir
