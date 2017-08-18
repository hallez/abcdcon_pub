#' ---
#'  title: ABCDCon Control Analysis &mdash; PS without top voxels
#'  author: Halle R. Dimsdale-Zucker
#'  output:
#'    html_document:
#'      toc: true
#'      number_sections: true
#'      theme: spacelab
#' ---

# this script is essentially `pattern_similarity_no_outlier_trials_load_data_btwn_runs.R` tweaked to read in files created by `control_analysis_PS_remove_top_voxels.m`

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
num_vox <- 5 # Set the number of voxels to pull; should also be set in `control_analysis_drop_voxels.R`

#' ## Setup other variables
subjects <- c(1:30)
exclude_subjects <- c(3:5,15:17)
subjects <- subjects[!is.element(subjects, exclude_subjects)]

rois <- c("CA1_body.nii", "CA2_3_DG_body.nii")

# initialize empty data frame
pattern_averages <- data.frame()
all_subj_patterns <- data.frame()
spatial_trials <- data.frame()
temporal_trials <- data.frame()

#' ### Flags
GRAPH_FLAG <- 0 # this will greatly slow down the script, but is nice for error-checking
MIXED_MODELS_FLAG <- 1
TALLY_CHECK_FLAG <- 1

#' # Loop and load in the data
#+ warning = FALSE
time_group_start <- Sys.time()
for(idir in c(1:length(roi_dirs))) {

  for(iroi in c(1:length(rois))){

    for(isub in subjects) {

      time_subj_start <- Sys.time()

      # define directories and filenames
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
      # this file is created by `control_analysis_PS_remove_top_voxels.m`
      cur_file <- file.path(
        cur_roi_dir,sprintf('br%s_pattern_corr_no_outlier_trials_%02d_truncated_voxels_all_runs.mat', cur_roi, num_vox))

      if(file.exists(cur_file)){

        # if want to graphically check after any step:
        # ggplot2::ggplot(reshape2::melt(subj_pattern_corr), ggplot2::aes(Var2, Var1, fill = value)) + ggplot2::geom_tile() + ggplot2::ggtitle(paste0(subj_number, " ", hemi_label, " ", cur_roi))

        # read in the subject's data for the current ROI
        subj_pattern_corr <- data.matrix(as.data.frame(R.matlab::readMat(cur_file)))

        # also read in the labels for the current trials
        subj_pattern_ids <- as.data.frame(R.matlab::readMat(paste0(
          cur_roi_dir,'br',cur_roi,'_pattern_mtx_ids_no_outlier_trials_all_runs.mat')))

        # TODO: figure out why this throws an error (even though output seems reasonable)
        subj_pattern_ids %>%
          tidyr::gather() -> subj_pattern_ids

        # capitalize on the power of R and give the pattern matrix meaningful
        # row and column names
        colnames(subj_pattern_corr) <- subj_pattern_ids$value
        rownames(subj_pattern_corr) <- subj_pattern_ids$value

        # notch out w/in run autocorrelation
        subj_pattern_corr[grep("Run[_]01",rownames(subj_pattern_corr)),
                          grep("Run[_]01",colnames(subj_pattern_corr))] <- NA
        subj_pattern_corr[grep("Run[_]02",rownames(subj_pattern_corr)),
                          grep("Run[_]02",colnames(subj_pattern_corr))] <- NA
        subj_pattern_corr[grep("Run[_]03",rownames(subj_pattern_corr)),
                          grep("Run[_]03",colnames(subj_pattern_corr))] <- NA
        subj_pattern_corr[grep("Run[_]04",rownames(subj_pattern_corr)),
                          grep("Run[_]04",colnames(subj_pattern_corr))] <- NA

        if(GRAPH_FLAG==2){
          # clear out temporary variables that have used previously
          tmp <- NULL
          tmp2 <- NULL

          # clean up labels
          subj_pattern_corr_for_graphing <- subj_pattern_corr

          # eliminate the room labels
          colnames(subj_pattern_corr_for_graphing) <- sub("Rm[(12)]","",colnames(subj_pattern_corr_for_graphing))
          rownames(subj_pattern_corr_for_graphing) <- sub("Rm[(12)]","",rownames(subj_pattern_corr_for_graphing))

          # reorder the data
          tmp2<-subj_pattern_corr_for_graphing[sort(rownames(subj_pattern_corr_for_graphing)),sort(colnames(subj_pattern_corr_for_graphing))]
          subj_pattern_corr_for_graphing <- tmp2

          # plot trial x trial correlation matrices
          print(ggplot2::ggplot(reshape2::melt(subj_pattern_corr_for_graphing), ggplot2::aes(Var2, Var1, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::ggtitle(paste0("Trial x Trial Correlations for ",
                           "subject:", subj_number,
                           "\nROI: ",cur_roi,
                           "\nData: ", fake_data_label)) +
            ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size=10,angle = 90, hjust = 1),
                  axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size=10)))
        } #GRAPH_FLAG == 2

        #' ## Extract trials of interest
        # clean up first
        # NB: will throw a warning the first time this is run b/c these variables won't exist
        remove(brown_brown_videos, gray_gray_videos,
               brown_house, gray_house,
               brown_rm1, brown_rm2, brown_rm1x2,
               gray_rm1, gray_rm2, gray_rm1x2,
               brown_vs_gray_house, brown1_gray2, brown2_gray1,
               brown_vid1, brown_vid2, brown_vid3, brown_vid4, brown_vid5, brown_vid6,
               brown_vid7, brown_vid8, brown_vid9, brown_vid10, brown_vid11, brown_vid12,
               brown_vid1_rm1, brown_vid2_rm1, brown_vid3_rm1, brown_vid4_rm1, brown_vid5_rm1, brown_vid6_rm1,
               brown_vid7_rm1, brown_vid8_rm1, brown_vid9_rm1, brown_vid10_rm1, brown_vid11_rm1, brown_vid12_rm1,
               brown_vid1_rm2, brown_vid2_rm2, brown_vid3_rm2, brown_vid4_rm2, brown_vid5_rm2, brown_vid6_rm2,
               brown_vid7_rm2, brown_vid8_rm2, brown_vid9_rm2, brown_vid10_rm2, brown_vid11_rm2, brown_vid12_rm2,
               gray_vid1, gray_vid2, gray_vid3, gray_vid4, gray_vid5, gray_vid6,
               gray_vid7, gray_vid8, gray_vid9, gray_vid10, gray_vid11, gray_vid12,
               gray_vid1_rm1, gray_vid2_rm1, gray_vid3_rm1, gray_vid4_rm1, gray_vid5_rm1, gray_vid6_rm1,
               gray_vid7_rm1, gray_vid8_rm1, gray_vid9_rm1, gray_vid10_rm1, gray_vid11_rm1, gray_vid12_rm1,
               gray_vid1_rm2, gray_vid2_rm2, gray_vid3_rm2, gray_vid4_rm2, gray_vid5_rm2, gray_vid6_rm2,
               gray_vid7_rm2, gray_vid8_rm2, gray_vid9_rm2, gray_vid10_rm2, gray_vid11_rm2, gray_vid12_rm2,
               brown_house_mtx, gray_house_mtx,
               brown_vs_gray_house_mtx,
               brown_vid1_mtx, brown_vid2_mtx, brown_vid3_mtx, brown_vid4_mtx, brown_vid5_mtx, brown_vid6_mtx,
               brown_vid7_mtx, brown_vid8_mtx, brown_vid9_mtx, brown_vid10_mtx, brown_vid11_mtx, brown_vid12_mtx,
               gray_vid1_mtx, gray_vid2_mtx, gray_vid3_mtx, gray_vid4_mtx, gray_vid5_mtx, gray_vid6_mtx,
               gray_vid7_mtx, gray_vid8_mtx, gray_vid9_mtx, gray_vid10_mtx, gray_vid11_mtx, gray_vid12_mtx,
               not_brown_vid1_mtx, not_brown_vid2_mtx, not_brown_vid3_mtx, not_brown_vid4_mtx, not_brown_vid5_mtx, not_brown_vid6_mtx,
               not_brown_vid7_mtx, not_brown_vid8_mtx, not_brown_vid9_mtx, not_brown_vid10_mtx, not_brown_vid11_mtx, not_brown_vid12_mtx,
               not_gray_vid1_mtx, not_gray_vid2_mtx, not_gray_vid3_mtx, not_gray_vid4_mtx, not_gray_vid5_mtx, not_gray_vid6_mtx,
               not_gray_vid7_mtx, not_gray_vid8_mtx, not_gray_vid9_mtx, not_gray_vid10_mtx, not_gray_vid11_mtx, not_gray_vid12_mtx)

        # check which objects are loaded
        objects()

        # sameHouse_sameRoom
        brown_rm1 <- subj_pattern_corr[grep("Brown[Rm1]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                       grep("Brown[Rm1]*[_]RHit[_]HC",colnames(subj_pattern_corr))]
        if(MIXED_MODELS_FLAG){
          brown_rm1_mtx <-
            brown_rm1 %>%
            as.data.frame() %>%
            dplyr::mutate(subj = subj_number,
                          roi = as.character(cur_roi),
                          hemi = hemi_label,
                          condition = "sameHouse_sameRoom",
                          row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

          if(TALLY_CHECK_FLAG){
            brown_rm1_mtx %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()
          }
        }

        brown_rm2 <- subj_pattern_corr[grep("Brown[Rm2]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                       grep("Brown[Rm2]*[_]RHit[_]HC",colnames(subj_pattern_corr))]
        if(MIXED_MODELS_FLAG){
          brown_rm2_mtx <-
            brown_rm2 %>%
            as.data.frame() %>%
            dplyr::mutate(subj = subj_number,
                          roi = as.character(cur_roi),
                          hemi = hemi_label,
                          condition = "sameHouse_sameRoom",
                          row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

          if(TALLY_CHECK_FLAG){
            brown_rm2_mtx %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()
          }
        }

        gray_rm1 <- subj_pattern_corr[grep("Gray[Rm1]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm1]*[_]RHit[_]HC",colnames(subj_pattern_corr))]
        if(MIXED_MODELS_FLAG){
          gray_rm1_mtx <-
            gray_rm1 %>%
            as.data.frame() %>%
            dplyr::mutate(subj = subj_number,
                          roi = as.character(cur_roi),
                          hemi = hemi_label,
                          condition = "sameHouse_sameRoom",
                          row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

          if(TALLY_CHECK_FLAG){
            gray_rm1_mtx %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()
          }
        }

        gray_rm2 <- subj_pattern_corr[grep("Gray[Rm2]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm2]*[_]RHit[_]HC",colnames(subj_pattern_corr))]
        if(MIXED_MODELS_FLAG){
          gray_rm2_mtx <-
            gray_rm2 %>%
            as.data.frame() %>%
            dplyr::mutate(subj = subj_number,
                          roi = as.character(cur_roi),
                          hemi = hemi_label,
                          condition = "sameHouse_sameRoom",
                          row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

          if(TALLY_CHECK_FLAG){
            gray_rm2_mtx %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()
          }
        }

        # sameHouse_differentRoom
        brown_rm1x2 <- subj_pattern_corr[grep("Brown[Rm1]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                         grep("Brown[Rm2]*[_]RHit[_]HC",colnames(subj_pattern_corr))]

        if(MIXED_MODELS_FLAG){
          brown_rm1x2_mtx <-
            brown_rm1x2 %>%
            as.data.frame() %>%
            dplyr::mutate(subj = subj_number,
                          roi = as.character(cur_roi),
                          hemi = hemi_label,
                          condition = "sameHouse_differentRoom",
                          row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

          if(TALLY_CHECK_FLAG){
            brown_rm1x2_mtx %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()

          }
        }

        gray_rm1x2 <- subj_pattern_corr[grep("Gray[Rm1]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                         grep("Gray[Rm2]*[_]RHit[_]HC",colnames(subj_pattern_corr))]

        if(MIXED_MODELS_FLAG){
          gray_rm1x2_mtx <-
            gray_rm1x2 %>%
            as.data.frame() %>%
            dplyr::mutate(subj = subj_number,
                          roi = as.character(cur_roi),
                          hemi = hemi_label,
                          condition = "sameHouse_differentRoom",
                          row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

          if(TALLY_CHECK_FLAG){
            gray_rm1x2_mtx %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()

          }
        }

        # sameHouse_eitherRoom
        brown_11_22_mtx <- dplyr::full_join(brown_rm1_mtx, brown_rm2_mtx, by = intersect(names(brown_rm1_mtx), names(brown_rm2_mtx)))
        brown_house_mtx <- dplyr::full_join(brown_11_22_mtx, brown_rm1x2_mtx, by = intersect(names(brown_11_22_mtx), names(brown_rm1x2_mtx))) %>%
          # squash prior values for `condition`
          dplyr::mutate(condition = "sameHouse_eitherRoom")

        if(TALLY_CHECK_FLAG){
          brown_house_mtx %>%
            tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
            tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
            dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
            dplyr::group_by(col_houseRoom, row_houseRoom) %>%
            dplyr::tally()
        }

        gray_11_22_mtx <- dplyr::full_join(gray_rm1_mtx, gray_rm2_mtx, by = intersect(names(gray_rm1_mtx), names(gray_rm2_mtx)))
        gray_house_mtx <- dplyr::full_join(gray_11_22_mtx, gray_rm1x2_mtx, by = intersect(names(gray_11_22_mtx), names(gray_rm1x2_mtx))) %>%
          # squash prior values for `condition`
          dplyr::mutate(condition = "sameHouse_eitherRoom")

        if(TALLY_CHECK_FLAG){
          gray_house_mtx %>%
            tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
            tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
            dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
            dplyr::group_by(col_houseRoom, row_houseRoom) %>%
            dplyr::tally()
        }

        # differentHouse_eitherRoom
        brown_vs_gray_house <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                                 grep("Brown[Rm12]*[_]RHit[_]HC",colnames(subj_pattern_corr))]

        if(MIXED_MODELS_FLAG==1) {
          brown_vs_gray_house_mtx <- brown_vs_gray_house %>%
            as.data.frame() %>%
            dplyr::mutate(subj = subj_number,
                          roi = as.character(cur_roi),
                          hemi = hemi_label,
                          condition = "differentHouse_eitherRoom",
                          row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

          if(TALLY_CHECK_FLAG){
            brown_vs_gray_house_mtx %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()
          }
        } #if MIXED_MODELS_FLAG

        # differentHouse_differentRoom
        brown1_gray2 <- subj_pattern_corr[grep("Brown[Rm1]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                          grep("Gray[Rm2]*[_]RHit[_]HC",colnames(subj_pattern_corr))]
        if(TALLY_CHECK_FLAG){
          brown1_gray2 %>%
            as.data.frame() %>%
            dplyr::mutate(row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name) %>%
            tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
            tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
            dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
            dplyr::group_by(col_houseRoom, row_houseRoom) %>%
            dplyr::tally()
        }

        brown2_gray1 <- subj_pattern_corr[grep("Brown[Rm2]*[_]RHit[_]HC",rownames(subj_pattern_corr)),
                                          grep("Gray[Rm1]*[_]RHit[_]HC",colnames(subj_pattern_corr))]
        if(TALLY_CHECK_FLAG){
          brown2_gray1 %>%
            as.data.frame() %>%
            dplyr::mutate(row_name = rownames(.)) %>%
            tidyr::gather(col_name, r, -row_name) %>%
            tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
            tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
            dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
            dplyr::group_by(col_houseRoom, row_houseRoom) %>%
            dplyr::tally()
        }

        # sameVideo, same room
        # for video1, need to make sure not
        # matching 10, 11, or 12 (hence $)
        # brown, room 1
        brown_vid1_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown1$",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown1$",colnames(subj_pattern_corr))]
        brown_vid2_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown2",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown2",colnames(subj_pattern_corr))]
        brown_vid3_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown3",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown3",colnames(subj_pattern_corr))]
        brown_vid4_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown4",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown4",colnames(subj_pattern_corr))]
        brown_vid5_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown5",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown5",colnames(subj_pattern_corr))]
        brown_vid6_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown6",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown6",colnames(subj_pattern_corr))]
        brown_vid7_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown7",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown7",colnames(subj_pattern_corr))]
        brown_vid8_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown8",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown8",colnames(subj_pattern_corr))]
        brown_vid9_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown9",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown9",colnames(subj_pattern_corr))]
        brown_vid10_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown10",rownames(subj_pattern_corr)),
                                             grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown10",colnames(subj_pattern_corr))]
        brown_vid11_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown11",rownames(subj_pattern_corr)),
                                             grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown11",colnames(subj_pattern_corr))]
        brown_vid12_rm1 <- subj_pattern_corr[grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown12",rownames(subj_pattern_corr)),
                                             grep("BrownRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Brown12",colnames(subj_pattern_corr))]

        # brown, room 2
        brown_vid1_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown1$",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown1$",colnames(subj_pattern_corr))]
        brown_vid2_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown2",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown2",colnames(subj_pattern_corr))]
        brown_vid3_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown3",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown3",colnames(subj_pattern_corr))]
        brown_vid4_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown4",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown4",colnames(subj_pattern_corr))]
        brown_vid5_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown5",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown5",colnames(subj_pattern_corr))]
        brown_vid6_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown6",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown6",colnames(subj_pattern_corr))]
        brown_vid7_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown7",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown7",colnames(subj_pattern_corr))]
        brown_vid8_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown8",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown8",colnames(subj_pattern_corr))]
        brown_vid9_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown9",rownames(subj_pattern_corr)),
                                            grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown9",colnames(subj_pattern_corr))]
        brown_vid10_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown10",rownames(subj_pattern_corr)),
                                             grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown10",colnames(subj_pattern_corr))]
        brown_vid11_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown11",rownames(subj_pattern_corr)),
                                             grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown11",colnames(subj_pattern_corr))]
        brown_vid12_rm2 <- subj_pattern_corr[grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown12",rownames(subj_pattern_corr)),
                                             grep("BrownRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Brown12",colnames(subj_pattern_corr))]

        # gray, room 1
        gray_vid1_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray1$",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray1$",colnames(subj_pattern_corr))]
        gray_vid2_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray2",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray2",colnames(subj_pattern_corr))]
        gray_vid3_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray3",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray3",colnames(subj_pattern_corr))]
        gray_vid4_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray4",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray4",colnames(subj_pattern_corr))]
        gray_vid5_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray5",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray5",colnames(subj_pattern_corr))]
        gray_vid6_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray6",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray6",colnames(subj_pattern_corr))]
        gray_vid7_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray7",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray7",colnames(subj_pattern_corr))]
        gray_vid8_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray8",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray8",colnames(subj_pattern_corr))]
        gray_vid9_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray9",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray9",colnames(subj_pattern_corr))]
        gray_vid10_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray10",rownames(subj_pattern_corr)),
                                             grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray10",colnames(subj_pattern_corr))]
        gray_vid11_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray11",rownames(subj_pattern_corr)),
                                             grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray11",colnames(subj_pattern_corr))]
        gray_vid12_rm1 <- subj_pattern_corr[grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray12",rownames(subj_pattern_corr)),
                                             grep("GrayRm[(1)][_]RHit[_]HCRC[_]EncVideo[_]Gray12",colnames(subj_pattern_corr))]

        # gray, room 2
        gray_vid1_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray1$",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray1$",colnames(subj_pattern_corr))]
        gray_vid2_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray2",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray2",colnames(subj_pattern_corr))]
        gray_vid3_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray3",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray3",colnames(subj_pattern_corr))]
        gray_vid4_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray4",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray4",colnames(subj_pattern_corr))]
        gray_vid5_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray5",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray5",colnames(subj_pattern_corr))]
        gray_vid6_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray6",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray6",colnames(subj_pattern_corr))]
        gray_vid7_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray7",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray7",colnames(subj_pattern_corr))]
        gray_vid8_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray8",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray8",colnames(subj_pattern_corr))]
        gray_vid9_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray9",rownames(subj_pattern_corr)),
                                            grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray9",colnames(subj_pattern_corr))]
        gray_vid10_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray10",rownames(subj_pattern_corr)),
                                             grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray10",colnames(subj_pattern_corr))]
        gray_vid11_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray11",rownames(subj_pattern_corr)),
                                             grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray11",colnames(subj_pattern_corr))]
        gray_vid12_rm2 <- subj_pattern_corr[grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray12",rownames(subj_pattern_corr)),
                                             grep("GrayRm[(2)][_]RHit[_]HCRC[_]EncVideo[_]Gray12",colnames(subj_pattern_corr))]

        if(TALLY_CHECK_FLAG){
          # this is commented out b/c would also need to handle the case where each df may not exist (eg, if a subject didn't see that video)
          #           brown_vid12_rm1 %>%
          #             as.data.frame() %>%
          #             dplyr::mutate(row_name = rownames(.)) %>%
          #             tidyr::gather(col_name, r, -row_name) %>%
          #             tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
          #             tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
          #             dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
          #             dplyr::group_by(col_houseRoom, col_vidID, row_houseRoom, row_vidID) %>%
          #             dplyr::tally()
        }

        # sameVideo
        # for video1, need to make sure not
        # matching 10, 11, or 12 (hence $)
        # enforce Rhits & house source hits (but NOT room hits)
        brown_vid1 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown1$",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown1$",colnames(subj_pattern_corr))]
        brown_vid2 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown2",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown2",colnames(subj_pattern_corr))]
        brown_vid3 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown3",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown3",colnames(subj_pattern_corr))]
        brown_vid4 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown4",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown4",colnames(subj_pattern_corr))]
        brown_vid5 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown5",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown5",colnames(subj_pattern_corr))]
        brown_vid6 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown6",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown6",colnames(subj_pattern_corr))]
        brown_vid7 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown7",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown7",colnames(subj_pattern_corr))]
        brown_vid8 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown8",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown8",colnames(subj_pattern_corr))]
        brown_vid9 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown9",rownames(subj_pattern_corr)),
                                        grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown9",colnames(subj_pattern_corr))]
        brown_vid10 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown10",rownames(subj_pattern_corr)),
                                         grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown10",colnames(subj_pattern_corr))]
        brown_vid11 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown11",rownames(subj_pattern_corr)),
                                         grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown11",colnames(subj_pattern_corr))]
        brown_vid12 <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown12",rownames(subj_pattern_corr)),
                                         grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown12",colnames(subj_pattern_corr))]

        gray_vid1 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray1$",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray1$",colnames(subj_pattern_corr))]
        gray_vid2 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray2",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray2",colnames(subj_pattern_corr))]
        gray_vid3 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray3",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray3",colnames(subj_pattern_corr))]
        gray_vid4 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray4",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray4",colnames(subj_pattern_corr))]
        gray_vid5 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray5",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray5",colnames(subj_pattern_corr))]
        gray_vid6 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray6",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray6",colnames(subj_pattern_corr))]
        gray_vid7 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray7",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray7",colnames(subj_pattern_corr))]
        gray_vid8 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray8",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray8",colnames(subj_pattern_corr))]
        gray_vid9 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray9",rownames(subj_pattern_corr)),
                                       grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray9",colnames(subj_pattern_corr))]
        gray_vid10 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray10",rownames(subj_pattern_corr)),
                                        grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray10",colnames(subj_pattern_corr))]
        gray_vid11 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray11",rownames(subj_pattern_corr)),
                                        grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray11",colnames(subj_pattern_corr))]
        gray_vid12 <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray12",rownames(subj_pattern_corr)),
                                        grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray12",colnames(subj_pattern_corr))]

        if(MIXED_MODELS_FLAG==1) {
          if(length(brown_vid1)>1){
            brown_vid1_mtx <- brown_vid1 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid1_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid2)>1){
            brown_vid2_mtx <- brown_vid2 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid2_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid3)>1){
            brown_vid3_mtx <- brown_vid3 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid3_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid4)>1){
            brown_vid4_mtx <- brown_vid4 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid4_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid5)>1){
            brown_vid5_mtx <- brown_vid5 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid5_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid6)>1){
            brown_vid6_mtx <- brown_vid6 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid6_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid7)>1){
            brown_vid7_mtx <- brown_vid7 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid7_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid8)>1){
            brown_vid8_mtx <- brown_vid8 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid8_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid9)>1){
            brown_vid9_mtx <- brown_vid9 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid9_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid10)>1){
            brown_vid10_mtx <- brown_vid10 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid10_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid11)>1){
            brown_vid11_mtx <- brown_vid11 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid11_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(brown_vid12)>1){
            brown_vid12_mtx <- brown_vid12 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              brown_vid12_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }


          if(length(gray_vid1)>1){
            gray_vid1_mtx <- gray_vid1 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid1_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid2)>1){
            gray_vid2_mtx <- gray_vid2 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid2_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid3)>1){
            gray_vid3_mtx <- gray_vid3 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid3_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid4)>1){
            gray_vid4_mtx <- gray_vid4 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid4_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid5)>1){
            gray_vid5_mtx <- gray_vid5 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid5_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid6)>1){
            gray_vid6_mtx <- gray_vid6 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid6_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid7)>1){
            gray_vid7_mtx <- gray_vid7 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid7_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid8)>1){
            gray_vid8_mtx <- gray_vid8 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid8_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid9)>1){
            gray_vid9_mtx <- gray_vid9 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid9_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid10)>1){
            gray_vid10_mtx <- gray_vid10 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid10_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid11)>1){
            gray_vid11_mtx <- gray_vid11 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid11_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }

          if(length(gray_vid12)>1){
            gray_vid12_mtx <- gray_vid12 %>%
              as.data.frame() %>%
              dplyr::mutate(subj = subj_number,
                            roi = as.character(cur_roi),
                            hemi = hemi_label,
                            condition = "sameVideo",
                            row_name = rownames(.)) %>%
              tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi)

            if(TALLY_CHECK_FLAG){
              gray_vid12_mtx %>%
                as.data.frame() %>%
                tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
                tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
                dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
                dplyr::group_by(col_vidID, row_vidID) %>%
                dplyr::tally()
            }
          }
        } #if MIXED_MODELS_FLAG

        # different video, same house
        brown_brown_videos <- subj_pattern_corr[grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown*",rownames(subj_pattern_corr)),
                                            grep("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown*",colnames(subj_pattern_corr))]

        brown_brown_videos_no_same <- NULL
        brown_brown_videos_no_same <- brown_brown_videos %>%
          as.data.frame() %>%
          dplyr::mutate(subj = subj_number,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        row_name = rownames(.)) %>%
          tidyr::gather(col_name, r, -row_name, -subj, -roi, -hemi) %>%
          # based on: http://stackoverflow.com/questions/32829358/dplyr-filter-with-sql-like-wildcard
          # http://stackoverflow.com/questions/4935479/how-to-combine-multiple-conditions-to-subset-a-data-frame-using-or
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown1$", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown1$", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown2", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown2", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown3", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown3", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown4", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown4", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown5", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown5", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown6", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown6", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown7", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown7", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown8", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown8", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown9", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown9", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown10", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown10", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown11", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown11", row_name))) %>%
          dplyr::filter(!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown12", col_name) | (!grepl("Brown[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Brown12", row_name))) %>%
          dplyr::mutate(condition = "differentVideo_sameHouse")

        head(brown_brown_videos_no_same)

        if(TALLY_CHECK_FLAG){
          brown_brown_videos_no_same_spread <- NULL
          brown_brown_videos_no_same_spread <-
            brown_brown_videos_no_same %>%
            tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
            tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
            dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel)

          # use this formatting so can see that there are zeros along the diagonal
          table(brown_brown_videos_no_same_spread$col_vidID, brown_brown_videos_no_same_spread$row_vidID)
        }


        gray_gray_videos <- subj_pattern_corr[grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray*",rownames(subj_pattern_corr)),
                                              grep("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray*",colnames(subj_pattern_corr))]

        gray_gray_videos_no_same <- gray_gray_videos %>%
          as.data.frame() %>%
          dplyr::mutate(subj = subj_number,
                        roi = as.character(cur_roi),
                        hemi = hemi_label,
                        condition = "gray_gray_video",
                        row_name = rownames(.)) %>%
          tidyr::gather(col_name, r, -row_name, -condition, -subj, -roi, -hemi) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray1$", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray1$", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray2", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray2", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray3", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray3", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray4", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray4", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray5", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray5", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray6", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray6", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray7", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray7", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray8", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray8", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray9", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray9", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray10", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray10", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray11", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray11", row_name))) %>%
          dplyr::filter(!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray12", col_name) | (!grepl("Gray[Rm12]*[_]RHit[_]HCR[CI]*[_]EncVideo[_]Gray12", row_name))) %>%
          dplyr::mutate(condition = "differentVideo_sameHouse")

        head(gray_gray_videos_no_same)

        if(TALLY_CHECK_FLAG){
          gray_gray_videos_no_same_spread <- NULL
          gray_gray_videos_no_same_spread <-
            gray_gray_videos_no_same %>%
            tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
            tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
            dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel)

          table(gray_gray_videos_no_same_spread$col_vidID, gray_gray_videos_no_same_spread$row_vidID)
        }


        # take averages - this is useful to save out if want to run an ANOVA instead of a mixed model
        if(dim(pattern_averages)[1]==0){
          pattern_averages[1,1]$subj <- subj_number
          pattern_averages[1,2]$roi <- as.character(cur_roi)
          pattern_averages[1,3]$hemi <- hemi_label
          pattern_averages[1,4]$same_house <- mean(c(brown_house_mtx$r, gray_house_mtx$r),
                                                   na.rm=TRUE)
          pattern_averages[1,5]$different_house <- mean(brown_vs_gray_house, na.rm=TRUE)
          pattern_averages[1,6]$sameHouse_sameRoom <- mean(c(brown_rm1, brown_rm2,
                                                             gray_rm1, gray_rm2),
                                                           na.rm=TRUE)
          pattern_averages[1,7]$sameHouse_differentRoom <- mean(c(brown_rm1x2,
                                                                  gray_rm1x2),
                                                                na.rm=TRUE)
          pattern_averages[1,8]$differentHouse_differentRoom <- mean(c(brown1_gray2,
                                                                       brown2_gray1),
                                                                     na.rm=TRUE)
          pattern_averages[1,9]$sameVideo <- mean(c(brown_vid1,brown_vid2,brown_vid3,brown_vid4,brown_vid5,brown_vid6,
                                                    brown_vid7,brown_vid8,brown_vid9,brown_vid10,brown_vid11,brown_vid12,
                                                    gray_vid1,gray_vid2,gray_vid3,gray_vid4,gray_vid5,gray_vid6,
                                                    gray_vid7,gray_vid8,gray_vid9,gray_vid10,gray_vid11,gray_vid12),
                                                  na.rm=TRUE)
          pattern_averages[1,10]$differentVideo_sameHouse <- mean(c(brown_brown_videos_no_same$r, gray_gray_videos_no_same$r), na.rm=TRUE)
          pattern_averages[1,11]$sameVideo_sameRoom <- mean(c(brown_vid1_rm1,brown_vid2_rm1,brown_vid3_rm1,brown_vid4_rm1,brown_vid5_rm1,brown_vid6_rm1,
                                                    brown_vid7_rm1,brown_vid8_rm1,brown_vid9_rm1,brown_vid10_rm1,brown_vid11_rm1,brown_vid12_rm1,
                                                    gray_vid1_rm1,gray_vid2_rm1,gray_vid3_rm1,gray_vid4_rm1,gray_vid5_rm1,gray_vid6_rm1,
                                                    gray_vid7_rm1,gray_vid8_rm1,gray_vid9_rm1,gray_vid10_rm1,gray_vid11_rm1,gray_vid12_rm1,
                                                    brown_vid1_rm2,brown_vid2_rm2,brown_vid3_rm2,brown_vid4_rm2,brown_vid5_rm2,brown_vid6_rm2,
                                                    brown_vid7_rm2,brown_vid8_rm2,brown_vid9_rm2,brown_vid10_rm2,brown_vid11_rm2,brown_vid12_rm2,
                                                    gray_vid1_rm2,gray_vid2_rm2,gray_vid3_rm2,gray_vid4_rm2,gray_vid5_rm2,gray_vid6_rm2,
                                                    gray_vid7_rm2,gray_vid8_rm2,gray_vid9_rm2,gray_vid10_rm2,gray_vid11_rm2,gray_vid12_rm2),
                                                  na.rm=TRUE)
        } else {
          pattern_averages[row_dimension,]$subj <- subj_number
          pattern_averages[row_dimension,]$roi <- as.character(cur_roi)
          pattern_averages[row_dimension,]$hemi <- hemi_label
          pattern_averages[row_dimension,]$same_house <- mean(c(brown_house_mtx$r, gray_house_mtx$r),
                                                              na.rm=TRUE)
          pattern_averages[row_dimension,]$different_house <- mean(brown_vs_gray_house, na.rm=TRUE)
          pattern_averages[row_dimension,]$sameHouse_sameRoom <- mean(c(brown_rm1, brown_rm2,
                                                                        gray_rm1, gray_rm2),
                                                                      na.rm=TRUE)
          pattern_averages[row_dimension,]$sameHouse_differentRoom <- mean(c(brown_rm1x2,
                                                                             gray_rm1x2),
                                                                           na.rm=TRUE)
          pattern_averages[row_dimension,]$differentHouse_differentRoom <- mean(c(brown1_gray2,
                                                                                  brown2_gray1),
                                                                                na.rm=TRUE)
          pattern_averages[row_dimension,]$sameVideo <- mean(c(brown_vid1,brown_vid2,brown_vid3,brown_vid4,brown_vid5,brown_vid6,
                                                               brown_vid7,brown_vid8,brown_vid9,brown_vid10,brown_vid11,brown_vid12,
                                                               gray_vid1,gray_vid2,gray_vid3,gray_vid4,gray_vid5,gray_vid6,
                                                               gray_vid7,gray_vid8,gray_vid9,gray_vid10,gray_vid11,gray_vid12),
                                                             na.rm=TRUE)
          pattern_averages[row_dimension,]$differentVideo_sameHouse <- mean(c(brown_brown_videos_no_same$r, gray_gray_videos_no_same$r), na.rm=TRUE)
          pattern_averages[row_dimension,]$sameVideo_sameRoom <- mean(c(brown_vid1_rm1,brown_vid2_rm1,brown_vid3_rm1,brown_vid4_rm1,brown_vid5_rm1,brown_vid6_rm1,
                                                              brown_vid7_rm1,brown_vid8_rm1,brown_vid9_rm1,brown_vid10_rm1,brown_vid11_rm1,brown_vid12_rm1,
                                                              gray_vid1_rm1,gray_vid2_rm1,gray_vid3_rm1,gray_vid4_rm1,gray_vid5_rm1,gray_vid6_rm1,
                                                              gray_vid7_rm1,gray_vid8_rm1,gray_vid9_rm1,gray_vid10_rm1,gray_vid11_rm1,gray_vid12_rm1,
                                                              brown_vid1_rm2,brown_vid2_rm2,brown_vid3_rm2,brown_vid4_rm2,brown_vid5_rm2,brown_vid6_rm2,
                                                              brown_vid7_rm2,brown_vid8_rm2,brown_vid9_rm2,brown_vid10_rm2,brown_vid11_rm2,brown_vid12_rm2,
                                                              gray_vid1_rm2,gray_vid2_rm2,gray_vid3_rm2,gray_vid4_rm2,gray_vid5_rm2,gray_vid6_rm2,
                                                              gray_vid7_rm2,gray_vid8_rm2,gray_vid9_rm2,gray_vid10_rm2,gray_vid11_rm2,gray_vid12_rm2),
                                                            na.rm=TRUE)
        }

        head(pattern_averages)
        row_dimension <- dim(pattern_averages)[1] + 1

        # clean up before moving onto next subject
        subj_pattern_ids <- NULL
        subj_pattern_corr <- NULL

        if(MIXED_MODELS_FLAG==1) {
          #' # Combine single-trial matrices for conditions of interest across subjects
          # NB: will get an error message about joining factors w/ different levels. this is okay.
          same_house_tmp <- NULL
          spatial_tmp <- NULL
          same_house_tmp <- dplyr::full_join(brown_house_mtx, gray_house_mtx, by = intersect(names(brown_house_mtx), names(gray_house_mtx)))
          # sanity check
          summary(same_house_tmp)
          dim(brown_house_mtx)
          dim(gray_house_mtx)
          # if we add these together, should get the dimensions of `same_house_tmp`
          dim(same_house_tmp)

          if(TALLY_CHECK_FLAG){
            same_house_tmp %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()
          }

          # continue combining
          spatial_tmp <- dplyr::full_join(same_house_tmp, brown_vs_gray_house_mtx, by = intersect(names(same_house_tmp), names(brown_vs_gray_house_mtx)))
          print(cat("dim(spatial_tmp)"))
          print(dim(spatial_tmp))
          # more sanity checking
          # and if add `dim(same_house_tmp)` and `dim(brown_vs_gray_house_mtx)` should get `dim(spatial_tmp)`
          dim(same_house_tmp)
          dim(brown_vs_gray_house_mtx)
          dim(spatial_tmp)

          if(TALLY_CHECK_FLAG){
            spatial_tmp %>%
              tidyr::separate(row_name, into = c("row_run", "row_run#", "row_trial", "row_trial#", "row_houseRoom", "row_objMem", "row_sourceMem", "row_vidLabel", "row_vidID")) %>%
              tidyr::separate(col_name, into = c("col_run", "col_run#", "col_trial", "col_trial#", "col_houseRoom", "col_objMem", "col_sourceMem", "col_vidLabel", "col_vidID")) %>%
              dplyr::select(-row_run, -row_trial, -row_vidLabel, -col_run, -col_trial, -col_vidLabel) %>%
              dplyr::group_by(col_houseRoom, row_houseRoom) %>%
              dplyr::tally()
          }

          # now we can eliminate NA rows
          # `dim(spatial_tmp_no_NA)` should be LESS than `dim(spatial_tmp)`
          spatial_tmp_no_NA <- spatial_tmp[!is.na(spatial_tmp$r),]
          dim(spatial_tmp_no_NA)
          print(cat("dim(spatial_tmp_no_NA)"))
          print(dim(spatial_tmp_no_NA))
          # clear out variables ASAP to free up space
          same_house_tmp <- NULL
          spatial_tmp <- NULL

          # even thought it's generally poor form to overwrite the same variable, doing it to deal w/ when one of the video_mtx variables doesn't exist (which is expected since not all trial conditions are RHits)
          # s021 is a good test case for when `brown_vid2_mtx` doesn't exist
          remove(same_video_tmp)

          if(exists("brown_vid1_mtx")){
            same_video_tmp <- brown_vid1_mtx
          }

          if(exists("same_video_tmp") && exists("brown_vid2_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid2_mtx, by = intersect(names(same_video_tmp), names(brown_vid2_mtx)))
          } else if (!exists("same_video_tmp") && exists("brown_vid2_mtx")) {
            same_video_tmp <- brown_vid2_mtx
          }

          if(exists("same_video_tmp") && exists("brown_vid3_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid3_mtx, by = intersect(names(same_video_tmp), names(brown_vid3_mtx)))
          } else if (!exists("same_video_tmp") && exists("brown_vid3_mtx")){
            same_video_tmp <- brown_vid3_mtx
          }

          if(exists("same_video_tmp") && exists("brown_vid4_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid4_mtx, by = intersect(names(same_video_tmp), names(brown_vid4_mtx)))
          } else if (!exists("same_video_tmp") && exists("brown_vid4_mtx")) {
            same_video_tmp <- brown_vid4_mtx
          }

          if(exists("brown_vid5_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid5_mtx, by = intersect(names(same_video_tmp), names(brown_vid5_mtx)))
          }

          if(exists("brown_vid6_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid6_mtx, by = intersect(names(same_video_tmp), names(brown_vid6_mtx)))
          }

          if(exists("brown_vid7_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid7_mtx, by = intersect(names(same_video_tmp), names(brown_vid7_mtx)))
          }

          if(exists("brown_vid8_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid8_mtx, by = intersect(names(same_video_tmp), names(brown_vid8_mtx)))
          }

          if(exists("brown_vid9_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid9_mtx, by = intersect(names(same_video_tmp), names(brown_vid9_mtx)))
          }

          if(exists("brown_vid10_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid10_mtx, by = intersect(names(same_video_tmp), names(brown_vid10_mtx)))
          }

          if(exists("brown_vid11_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid11_mtx, by = intersect(names(same_video_tmp), names(brown_vid11_mtx)))
          }

          if(exists("brown_vid12_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, brown_vid12_mtx, by = intersect(names(same_video_tmp), names(brown_vid12_mtx)))
          }

          if(exists("gray_vid1_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid1_mtx, by = intersect(names(same_video_tmp), names(gray_vid1_mtx)))
          }

          if(exists("gray_vid2_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid2_mtx, by = intersect(names(same_video_tmp), names(gray_vid2_mtx)))
          }

          if(exists("gray_vid3_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid3_mtx, by = intersect(names(same_video_tmp), names(gray_vid3_mtx)))
          }

          if(exists("gray_vid4_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid4_mtx, by = intersect(names(same_video_tmp), names(gray_vid4_mtx)))
          }

          if(exists("gray_vid5_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid5_mtx, by = intersect(names(same_video_tmp), names(gray_vid5_mtx)))
          }

          if(exists("gray_vid6_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid6_mtx, by = intersect(names(same_video_tmp), names(gray_vid6_mtx)))
          }

          if(exists("gray_vid7_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid7_mtx, by = intersect(names(same_video_tmp), names(gray_vid7_mtx)))
          }

          if(exists("gray_vid8_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid8_mtx, by = intersect(names(same_video_tmp), names(gray_vid8_mtx)))
          }

          if(exists("gray_vid9_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid9_mtx, by = intersect(names(same_video_tmp), names(gray_vid9_mtx)))
          }

          if(exists("gray_vid10_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid10_mtx, by = intersect(names(same_video_tmp), names(gray_vid10_mtx)))
          }

          if(exists("gray_vid11_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid11_mtx, by = intersect(names(same_video_tmp), names(gray_vid11_mtx)))
          }

          if(exists("gray_vid12_mtx")){
            same_video_tmp <- dplyr::full_join(same_video_tmp, gray_vid12_mtx, by = intersect(names(same_video_tmp), names(gray_vid12_mtx)))
          }

          # remove NA values to reduce the size a bit
          # but first, sanity check the size
          dim(same_video_tmp)
          print(cat("dim(same_video_tmp)"))
          print(dim(same_video_tmp))
          same_video_tmp_no_NA <- same_video_tmp[!is.na(same_video_tmp$r),]
          # and check again after removing NA values
          dim(same_video_tmp_no_NA)
          print(cat("dim(same_video_tmp_no_NA)"))
          print(dim(same_video_tmp_no_NA))
          same_video_tmp <- NULL

          diff_video_tmp <- NULL
          diff_video_tmp <- dplyr::full_join(brown_brown_videos_no_same, gray_gray_videos_no_same, by = intersect(names(brown_brown_videos_no_same), names(gray_gray_videos_no_same)))

          # remove NA values to reduce the size a bit
          # but first, sanity check the size
          # should have more rows than `same_video_tmp` because more likely that two trials will NOT be in the same video
          dim(diff_video_tmp)
          print(cat("dim(diff_video_tmp)"))
          print(dim(diff_video_tmp))
          diff_video_tmp_no_NA <- diff_video_tmp[!is.na(diff_video_tmp$r),]
          # and check again after removing NA values
          dim(diff_video_tmp_no_NA)
          print(cat("dim(diff_video_tmp_no_NA)"))
          print(dim(diff_video_tmp_no_NA))
          diff_video_tmp <- NULL

          temporal_tmp <- NULL
          temporal_tmp <- dplyr::full_join(same_video_tmp_no_NA, diff_video_tmp_no_NA, by = intersect(names(same_video_tmp_no_NA), names(diff_video_tmp_no_NA)))
          # sanity check ourselves
          dim(temporal_tmp)
          print(cat("dim(temporal_tmp)"))
          print(dim(temporal_tmp))

          # write out individual subject files b/c they're huge
          save(spatial_tmp_no_NA,file=paste0(cur_roi_dir,
                                             paste0(subj_number, "_",hemi_label, "_", cur_roi,"_spatial_trial_patterns.RData")))
          write.csv(spatial_tmp_no_NA,file=paste0(cur_roi_dir,
                                                  paste0(subj_number, "_",hemi_label, "_", cur_roi,"_spatial_trial_patterns.csv")))
          save(temporal_tmp,file=paste0(cur_roi_dir,
                                        paste0(subj_number, "_",hemi_label, "_", cur_roi,"_temporal_trial_patterns.RData")))
          write.csv(temporal_tmp,file=paste0(cur_roi_dir,
                                             paste0(subj_number, "_",hemi_label, "_", cur_roi,"_temporal_trial_patterns.csv")))

          if(dim(spatial_trials)[1]==0 && dim(temporal_trials)[1]==0) {
            # deal w/ what happens for the first iteration
            spatial_trials <- spatial_tmp_no_NA
            temporal_trials <- temporal_tmp

            # clean up
            spatial_tmp_no_NA <- NULL
            temporal_tmp <- NULL
          } else {
            spatial_trials <- dplyr::full_join(spatial_trials, spatial_tmp_no_NA, by = intersect(names(spatial_tmp_no_NA), names(spatial_trials)))
            temporal_trials <- dplyr::full_join(temporal_trials, temporal_tmp, by = intersect(names(temporal_tmp), names(temporal_trials)))

            # clean up
            spatial_tmp_no_NA <- NULL
            temporal_tmp <- NULL
          } #if(isub == subjects[1])

          # let's do some quick sanity checking -- the number of rows should increase on each iteration
          print(cat("Let's check dimensions of `spatial_trials` and `temporal_trials` before iterating through next loop. "))
          print(dim(spatial_trials))
          print(dim(temporal_trials))
        } # if MIXED_MODELS_FLAG

        print(sprintf("It took %.2f seconds for %s %s %s.", Sys.time() - time_subj_start, subj_number, hemi_label, cur_roi))

      } else {
        print(sprintf('Pattern correlation file %s does not exist. Skipping.', cur_file))
        next
      }#if cur_file exists

    } #isub

  } #iroi

} #idir

# Print how long it took to go through all of the loops
sprintf("It took %.2f seconds to run all of this.", Sys.time() - time_group_start)

#' # Save out group data files
#' ## Take a quick peak at the data first
head(pattern_averages)
if(MIXED_MODELS_FLAG==1) {
  print(head(spatial_trials))
  print(head(temporal_trials))
}

#' # Print out included subjects
unique(pattern_averages$subj)

#' ## Actually save it out
# pattern averages
save(pattern_averages,file=paste0(analyzed_mri_dir, halle::ensure_trailing_slash("multivariate_sanityCheck"), sprintf('group_pattern_averages_no_outlier_trials_%02d_truncated_voxels_btwn_runs.RData', num_vox)))
write.csv(pattern_averages,file=paste0(analyzed_mri_dir, halle::ensure_trailing_slash("multivariate_sanityCheck"), sprintf('group_pattern_averages_no_outlier_trials_%02d_truncated_voxels_btwn_runs.csv', num_vox)))

# r values for individual trials and individual subjects
if(MIXED_MODELS_FLAG==1){
  save(spatial_trials,file=paste0(analyzed_mri_dir, halle::ensure_trailing_slash("multivariate_sanityCheck"), sprintf('group_spatial_trial_patterns_%02d_truncated_voxels.RData', num_vox)))
  write.csv(spatial_trials,file=paste0(analyzed_mri_dir, halle::ensure_trailing_slash("multivariate_sanityCheck"), sprintf('group_spatial_trial_patterns_%02d_truncated_voxels.csv', num_vox)))
  save(temporal_trials,file=paste0(analyzed_mri_dir, halle::ensure_trailing_slash("multivariate_sanityCheck"), sprintf('group_temporal_trial_patterns_%02d_truncated_voxels.RData', num_vox)))
  write.csv(temporal_trials,file=paste0(analyzed_mri_dir, halle::ensure_trailing_slash("multivariate_sanityCheck"), sprintf('group_temporal_trial_patterns_%02d_truncated_voxels.csv', num_vox)))
}

