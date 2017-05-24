#' ---
#'  title: ABCDCon Figure out trial labels
#'  author: Halle R. Dimsdale-Zucker
#'  output:
#'    html_document:
#'      toc: true
#'      number_sections: true
#'      theme: spacelab
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

#' ## Flags
SAVE_FLAG <- 1

#' ## Set paths as variables
analyzed_mri_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_mri))
raw_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$raw_behavioral))
analyzed_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_behavioral))

#' ## Setup other variables
subjects <- c(1:30)
exclude_subjects <- c(3:5,16,17)
subjects <- subjects[!is.element(subjects, exclude_subjects)]
num_runs <- 4 #TODO: don't hard code

#' # Loop and load in the data

for(isub in subjects) {

  # overwrite variables that may have changed between subjects
  run_ids <- c(1,2,3,4)

  # define directories and filenames
  subj_number <- halle::format_three_digit_subject_id(isub)
  cur_subj_dir <- paste0(analyzed_mri_dir, halle::ensure_trailing_slash(subj_number))
  cur_subj_multivariate_dir <- paste0(analyzed_mri_dir, halle::ensure_trailing_slash("multivariate_sanityCheck"))

  # overwrite values for subjects w/ different numbers of runs
  if(subj_number=="s003"){
    run_ids <- c(1,0,0,4)
  } else if (subj_number=="s015"){
    run_ids <- c(1,2,0,4)
  }

  # initialize an empty variable where the trial IDs for all runs will go
  all_runs <- data.frame()

  for(runidx in c(1:num_runs)){

    irun <- run_ids[runidx]

    # test to make sure the current file exists
    # if not, skip and continue on
    subj_pattern_ids_fname <- NULL
    subj_pattern_ids_fname <- paste0(cur_subj_multivariate_dir,'pattern_mtx_ids_run_',sprintf("%02d",irun),'.mat')

    if(file.exists(subj_pattern_ids_fname)){

      # read in the trial labels
      subj_pattern_ids <- NULL
      subj_pattern_ids <- as.data.frame(R.matlab::readMat(subj_pattern_ids_fname))

      # clean up the labels
      subj_pattern_ids_tidy <- NULL

      subj_pattern_ids_tidy <-
        subj_pattern_ids %>%
        # this will throw a warning about attributes being non-identical
        # as far as I can tell, it warns about dropping data but doesn't actually
        # for potential solution, see:
        # http://stackoverflow.com/questions/28972386/retain-attributes-when-using-gather-from-tidyr-attributes-are-not-identical
        tidyr::gather() %>%
        # nicely format Run##_Trial### for all rows (becuase this is common to each trial)
        dplyr::mutate(tidy_lbls = stringr::str_replace(value,"Run[_](.*?)[_]Trial[_](.*?)", "Run\\1_Trial\\2")) %>%
        # put a placeholder RoomID for hits w/o room labels
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "Brown_", "Brown_RoomID_")) %>%
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "Gray_", "Gray_RoomID_")) %>%
        # now, split out RoomID where it exists
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "BrownRm(.*?)_", "Brown_Rm\\1_")) %>%
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "GrayRm(.*?)_", "Gray_Rm\\1_")) %>%
        # now deal w/ all non-"R" response types
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls,"CR_EncVideo_(.*?)", "HouseLbl_RoomID_CR_LocationResp_EncVideo\\1")) %>%
        # this will include FA.R, FA.F, etc.
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "FA_(.*?)_EncVideo_(.*?)", "HouseLbl_RoomID_FA.\\1_LocationResp_EncVideo\\2")) %>%
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "Miss_EncVideo_(.*?)", "HouseLbl_RoomID_Miss_LocationResp_EncVideo\\1")) %>%
        # this only deals w/ trials that are excluded based on memory (or at least behavioral) performance,
        # trials identified to be excluded based on single trial modelling are NOT marked here
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "ExcludeTrial_EncVideo_(.*?)", "HouseLbl_RoomID_ExcludeTrial_LocationResp_EncVideo\\1")) %>%
        # tidy up the EncVideo label
        dplyr::mutate(tidy_lbls = stringr::str_replace_all(tidy_lbls, "EncVideo_(.*?)", "EncVideo\\1")) %>%
        # now that each trial label has 7 parts, just split them up
        tidyr::separate(tidy_lbls, into = c("run_id", "trial_id", "house_id", "room_id", "mem_resp_scored", "location_resp_scored", "enc_video_id"), sep = "_") %>%
        # change run_id and trial_id to numeric
        # these are eas(ier) because just have to remove verbal label and convert to numeric
        dplyr::mutate(run_id_numeric = stringr::str_replace_all(run_id,"Run(.*?)", "\\1")) %>%
        dplyr::mutate(run_id_numeric = as.numeric(run_id_numeric)) %>%
        dplyr::mutate(trial_id_numeric = stringr::str_replace_all(trial_id, "Trial(.*?)", "\\1")) %>%
        dplyr::mutate(trial_id_numeric = as.numeric(trial_id_numeric)) %>%
        # convert house_id to numeric
        # 0 = no house label, 1 = brown house, 2 = gray house (it's just alphabetical)
        dplyr::mutate(house_id_numeric = stringr::str_replace_all(house_id, "HouseLbl", "0")) %>%
        # remember once start replacing values, need to overwrite w/in this new variable
        dplyr::mutate(house_id_numeric = stringr::str_replace_all(house_id_numeric, "Brown", "1")) %>%
        dplyr::mutate(house_id_numeric = stringr::str_replace_all(house_id_numeric, "Gray", "2")) %>%
        dplyr::mutate(house_id_numeric = as.numeric(house_id_numeric)) %>%
        # convert room_id to numeric
        # 0 = no room label, 1 = room1, 2 = room2
        dplyr::mutate(room_id_numeric = stringr::str_replace_all(room_id, "RoomID", "0")) %>%
        dplyr::mutate(room_id_numeric = stringr::str_replace_all(room_id_numeric, "Rm1", "1")) %>%
        dplyr::mutate(room_id_numeric = stringr::str_replace_all(room_id_numeric, "Rm2", "2")) %>%
        dplyr::mutate(room_id_numeric = as.numeric(room_id_numeric)) %>%
        # convert mem_resp_scored to numeric
        # 1 = RHit, 2 = FHit, 3 = CR, 4 = Miss, FA.R = 5, FA.F = 6, 99 = ExcludeTrial
        dplyr::mutate(mem_resp_scored_numeric = stringr::str_replace_all(mem_resp_scored, "RHit", "1")) %>%
        dplyr::mutate(mem_resp_scored_numeric = stringr::str_replace_all(mem_resp_scored_numeric, "FHit", "2")) %>%
        dplyr::mutate(mem_resp_scored_numeric = stringr::str_replace_all(mem_resp_scored_numeric, "CR", "3")) %>%
        dplyr::mutate(mem_resp_scored_numeric = stringr::str_replace_all(mem_resp_scored_numeric, "Miss", "4")) %>%
        dplyr::mutate(mem_resp_scored_numeric = stringr::str_replace_all(mem_resp_scored_numeric, "FA.R", "5")) %>%
        dplyr::mutate(mem_resp_scored_numeric = stringr::str_replace_all(mem_resp_scored_numeric, "FA.F", "6")) %>%
        dplyr::mutate(mem_resp_scored_numeric = stringr::str_replace_all(mem_resp_scored_numeric, "ExcludeTrial", "99")) %>%
        dplyr::mutate(mem_resp_scored_numeric = as.numeric(mem_resp_scored_numeric)) %>%
        # convert location_resp_scored to numeric
        # 0 = no resp ("LocationResp"), 1 = HCRC (house correct, room correct), 2 = HCRI (house correct, room incorrect), 3 = HI (house incorrect)
        dplyr::mutate(location_resp_scored_numeric = stringr::str_replace_all(location_resp_scored, "LocationResp", "0")) %>%
        dplyr::mutate(location_resp_scored_numeric = stringr::str_replace_all(location_resp_scored_numeric, "HCRC", "1")) %>%
        dplyr::mutate(location_resp_scored_numeric = stringr::str_replace_all(location_resp_scored_numeric, "HCRI", "2")) %>%
        dplyr::mutate(location_resp_scored_numeric = stringr::str_replace_all(location_resp_scored_numeric, "HI", "3")) %>%
        dplyr::mutate(location_resp_scored_numeric = as.numeric(location_resp_scored_numeric)) %>%
        # convert enc_video_id to numeric
        # the first digit indicates house ID (1 = brown, 2 = gray)
        # the next 2 digits indicate video ID
        # ie, brown1 would be 101, gray1 would be 201
        # brown10 would be 110, gray12 would be 212
        # lure = 0
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id, "EncVideoLure", "0")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown2", "102")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown3", "103")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown4", "104")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown5", "105")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown6", "106")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown7", "107")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown8", "108")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown9", "109")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown10", "110")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown11", "111")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown12", "112")) %>%
        # if try to match on Brown1 before have done videos 10, 11, and 12 will end up with weird 4-digit values
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoBrown1", "101")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray2", "202")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray3", "203")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray4", "204")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray5", "205")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray6", "206")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray7", "207")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray8", "208")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray9", "209")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray10", "210")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray11", "211")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray12", "212")) %>%
        dplyr::mutate(enc_video_id_numeric = stringr::str_replace_all(enc_video_id_numeric, "EncVideoGray1", "201")) %>%
        dplyr::mutate(enc_video_id_numeric = as.numeric(enc_video_id_numeric)) %>%
        # select only relevant columns
        dplyr::select(-key, -value)

      # add a column that just indicates if a trial is excluded
      # first, mark trials that get excluded for behavioral performance
      # this is helpful in checking if number of rhits (and other relevant trials)
      # match w/ behaviorally scored data
      subj_pattern_ids_tidy <-
        subj_pattern_ids_tidy %>%
        dplyr::mutate(excluded_trials_numeric = ifelse(mem_resp_scored_numeric==99,1,0))

      # mark trials identified to be excluded from `beta_timeseries_id_outliers.m`
      cur_excluded_trials_fname <- paste0(cur_subj_multivariate_dir, halle::ensure_trailing_slash(subj_number), subj_number,"_bad_timepoints_run", sprintf("%02d",irun), ".dat")
      if(file.exists(cur_excluded_trials_fname)){
        cur_excluded_trials <- NULL
        cur_excluded_trials <- read.table(cur_excluded_trials_fname, header = FALSE, sep = ",")

        cur_excluded_trials_tidy <- NULL
        cur_excluded_trials_tidy <-
          cur_excluded_trials %>%
          tidyr::gather(col_name,excluded_trial_num)

        # overwrite `mem_resp_scored_numeric` to 99 (ie, the code for excluded trials)
        # based on: http://stackoverflow.com/questions/5824173/replace-a-value-in-a-data-frame-based-on-a-conditional-if-statement-in-r
        subj_pattern_ids_tidy$mem_resp_scored_numeric[subj_pattern_ids_tidy$trial_id_numeric %in% cur_excluded_trials_tidy$excluded_trial_num] <- 99
      } #if(file.exists(cur_excluded_trials_fname

      # now modify this column for trials that get excluded based on aberrant beta
      # so 1 = exclude for behavioral reasons
      # and 2 = exclude as aberrant beta
      # this is helpful when verifying if pulling correct trials
      subj_pattern_ids_tidy <-
        subj_pattern_ids_tidy %>%
        dplyr::mutate(excluded_trials_numeric = ifelse(mem_resp_scored_numeric==99 & excluded_trials_numeric==0,2,
                                                       ifelse(mem_resp_scored_numeric==99,1,0)))

      # put these data in a variable for all runs
      if(dim(all_runs)[1]==0){
        all_runs <- subj_pattern_ids_tidy
      } else {
        # merge data here
        all_runs <- dplyr::full_join(all_runs, subj_pattern_ids_tidy, by = intersect(names(all_runs), names(subj_pattern_ids_tidy)))
      }

    } #if cur_file exists

    # even if cleared earlier, can never be too safe!
    subj_pattern_ids <- NULL
    subj_pattern_ids_tidy <- NULL

  } #irun

  # eliminate any non-numeric only columns
  all_runs_for_matlab <- NULL
  all_runs_for_matlab <-
    all_runs %>%
    dplyr::select(contains("_numeric")) %>%
    # ensure trials are in order w/in a run
    # this is very, very important
    # b/c later will be indexing according to order of betas
    # (which are in ascending order of trials)
    dplyr::arrange(run_id_numeric,trial_id_numeric)

  # take a quick peak at the data first
  head(all_runs_for_matlab)

  # actually save it out
  if(SAVE_FLAG){
    fname_out <- paste0(cur_subj_multivariate_dir,halle::ensure_trailing_slash(subj_number),subj_number,'_trial_pattern_ids_all_runs.mat')
    R.matlab::writeMat(fname_out, ids = all_runs_for_matlab)
  }

  # clear variables before moving onto next subject
  # again, even if cleared earlier, better to be safe than sorry
  all_runs <- NULL
  all_runs_for_matlab <- NULL
} #isub
