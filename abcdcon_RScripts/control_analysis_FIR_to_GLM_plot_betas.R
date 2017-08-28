#' ---
#' title: ABCDCon FIR Control Analysis
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
graph_fpath_out <- paste0(dropbox_graph_fpath_out, "FIR-to-GLM")
dir.create(graph_fpath_out) # will throw an error if this already exists

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
ROIs <- c("brCA1_body", "brCA2_3_DG_body", "brwhole_hippo")
model_name <- "byMemoryFHitMiss"
all_betas <- data.frame()

#' # Load in data across hemis, ROIs, and subjects
for(ihemi in 1:length(hemis)){
  cur_hemi <- hemis[ihemi]

  for(iroi in 1:length(ROIs)){
    cur_roi <- ROIs[iroi]

    for(isubj in 1:length(subj_formatted)){
      cur_subj <- subj_formatted[isubj]

      cur_betas_fpath <- file.path(analyzed_mri_dir, cur_subj, 'ROIs', cur_hemi, sprintf('%s_%s_FIR_beta_nan_means_all_runs.mat', cur_roi, model_name))
      cur_ids_fpath <- file.path(analyzed_mri_dir, cur_subj, 'ROIs', cur_hemi, sprintf('%s_%s_FIR_beta_ids_all_runs.mat', cur_roi, model_name))
      if(file.exists(cur_betas_fpath)){
        cur_betas <- as.data.frame(R.matlab::readMat(cur_betas_fpath))
        cur_ids <- as.data.frame(R.matlab::readMat(cur_ids_fpath))

        # --- Label `cur_betas` w/ `cur_ids` ---
        # there's probably a better dplyr way to do this, but works (and is the same as in `pattern_similarity_no_outlier_trials_load_data_btwn_runs`)
        # will get the following warning - `attributes are not identical across measure variables; they will be dropped `
        cur_ids_fmt <- cur_ids %>%
          tidyr::gather()

        cur_betas_w_ids <- cur_betas
        colnames(cur_betas_w_ids) <- cur_ids_fmt$value

        # --- Flip dataframe around so have rows instead of cols ---
        cur_betas_tidy <- cur_betas_w_ids %>%
          tidyr::gather(key = regressor_type, value = mean_beta_val) %>%
          dplyr::mutate(subj_id = cur_subj,
                        roi = cur_roi,
                        hemi = cur_hemi)

        # --- Put into a group data frame ---
        if(dim(all_betas)[1]==0){
          all_betas <- cur_betas_tidy
        } else {
          all_betas <- dplyr::full_join(all_betas, cur_betas_tidy, by = intersect(names(all_betas), names(cur_betas_tidy)))
        } #if dim(all_betas

      } else {
        print(sprintf("Current mean betas file %s does not exist. Continuing on.", cur_betas_fpath))
        next
      } #if(file.exists
    }# for isubj
  } #for iroi
} #for ihemi

#' # Tidy up group betas dataframe
all_betas_tidy <- all_betas %>%
  tidyr::separate(regressor_type, into = c("regressor_name", "beta_ID", "run_number"), sep = "_") %>%
  dplyr::mutate(beta_number_bare = sub("beta", "", beta_ID)) %>%
  dplyr::mutate(beta_number = as.numeric(beta_number_bare)) %>%
  dplyr::mutate(beta_fact = as.factor(beta_number)) %>%
  # filter out regressors for motion, etc.
  dplyr::filter(regressor_name %in% c("RHit", "FHitsANDMisses", "FA", "CR")) %>% #NB: currently won't work for both FIR models (FHit, Miss vs FHitsANDMisses)
  # ensure that trials are in order (b/c sequencing betas assumes this)
  dplyr::arrange(subj_id, hemi, roi, run_number, beta_ID) %>%
  # based on https://stackoverflow.com/questions/30793033/r-add-columns-indicating-start-and-end-for-a-sequence-within-columns (see setup in question post)
  dplyr::group_by(subj_id, hemi, roi, regressor_name, run_number) %>%
  dplyr::mutate(beta_seq = seq_along(regressor_name)) %>%
  dplyr::mutate(beta_seq_fact = as.factor(beta_seq)) %>%
  dplyr::ungroup() %>%
  # create better ROI labels for plots
  dplyr::mutate(roi_lbl = ifelse(roi == "brCA1_body", "CA1",
                                 ifelse(roi == "brCA2_3_DG_body", "CA23DG",
                                        ifelse(roi == "brwhole_hippo", "HC", roi))))

#' ## Peek at tidied data
head(all_betas_tidy)
unique(all_betas_tidy$regressor_name)

#' # Print what summarized values should be
all_betas_tidy %>%
  dplyr::filter(regressor_name %in% c("RHit", "FHitsANDMisses")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq, regressor_name) %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val)

#' # Plot
#' ### both hemi
all_betas_tidy %>%
  dplyr::filter(regressor_name == "RHit") %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq) %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = roi_lbl)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = min_val, ymax = max_val)) +
  ggplot2::facet_grid(.~hemi) +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHits_geomline_both-hemi.pdf"),
                  width=8, height=6)
}

#' ### both hemi, shaded error bars
all_betas_tidy %>%
  dplyr::filter(regressor_name == "RHit") %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq) %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = roi_lbl)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = min_val, ymax = max_val, color = roi_lbl, fill = roi_lbl), alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(.~hemi) +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHits_geomline_both-hemi_shaded-errorbars.pdf"),
                  width=8, height=6)
}

#' ### both hemi, no errorbars
all_betas_tidy %>%
  dplyr::filter(regressor_name == "RHit") %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq) %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = roi_lbl)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_grid(.~hemi) +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHits_geomline_both-hemi_no-errorbars.pdf"),
                  width=8, height=6)
}

#' ### just left hemi
all_betas_tidy %>%
  dplyr::filter(regressor_name == "RHit") %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq) %>%
  dplyr::filter(hemi == "ashs_left") %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = roi_lbl)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = min_val, ymax = max_val))  +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHits_geomline_left-hemi.pdf"),
                  width=8, height=6)
}

#' ### just left hemi, shaded error bars
all_betas_tidy %>%
  dplyr::filter(regressor_name == "RHit") %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq) %>%
  dplyr::filter(hemi == "ashs_left") %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = roi_lbl)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = min_val, ymax = max_val, color = roi_lbl, fill = roi_lbl), alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHits_geomline_left-hemi_shaded-errorbars.pdf"),
                  width=8, height=6)
}

#' ### just left hemi, facet by roi_lbl
all_betas_tidy %>%
  dplyr::filter(regressor_name == "RHit") %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq) %>%
  dplyr::filter(hemi == "ashs_left") %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = roi_lbl)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = min_val, ymax = max_val, color = regressor_name, fill = regressor_name), alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(.~roi_lbl)  +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHits_geomline_left-hemi_facet-roi_shaded-errorbars.pdf"),
                  width=8, height=6)
}

#' ### R, FHitANDMiss: both hemi, shaded error bars
all_betas_tidy %>%
  dplyr::filter(regressor_name %in% c("RHit", "FHitsANDMisses")) %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::group_by(hemi, roi_lbl, beta_seq, regressor_name) %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = regressor_name)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = min_val, ymax = max_val, color = regressor_name, fill = regressor_name), alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(roi_lbl~hemi) +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHitsVSFhitsMiss_geomline_both-hemi_shaded-errorbars.pdf"),
                  width=8, height=6)
}

#' ### R, FHitANDMiss: left hemi, shaded error bars
all_betas_tidy %>%
  dplyr::filter(regressor_name %in% c("RHit", "FHitsANDMisses")) %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG")) %>%
  dplyr::filter(hemi == "ashs_left") %>%
  dplyr::group_by(roi_lbl, beta_seq, regressor_name) %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = regressor_name)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = min_val, ymax = max_val, color = regressor_name, fill = regressor_name), alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(roi_lbl~.) +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR


if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHitsVSFhitsMiss_geomline_left-hemi_shaded-errorbars.pdf"),
                  width=7, height=6)
}

#' ### R, FHitANDMiss: left hemi, shaded error bars, CA1/CA23DG (body) and HC
all_betas_tidy %>%
  dplyr::filter(regressor_name %in% c("RHit", "FHitsANDMisses")) %>%
  dplyr::filter(roi_lbl %in% c("CA1", "CA23DG", "HC")) %>%
  dplyr::filter(hemi == "ashs_left") %>%
  dplyr::group_by(roi_lbl, beta_seq, regressor_name) %>%
  dplyr::summarise(gmean_beta_val = mean(mean_beta_val, na.rm = TRUE),
                   num_obs = length(mean_beta_val), # should be 8 if the subject had one of each trial type per run
                   sem_beta_val = sd(mean_beta_val, na.rm = T) / sqrt(num_obs),
                   min_val = gmean_beta_val - sem_beta_val,
                   max_val = gmean_beta_val + sem_beta_val) %>%
  ggplot2::ggplot(ggplot2::aes(x = beta_seq, y = gmean_beta_val, color = regressor_name)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = min_val, ymax = max_val, color = regressor_name, fill = regressor_name), alpha = 0.2) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(roi_lbl~.) +
  ggplot2::ylab("mean beta value") +
  ggplot2::xlab("FIR timepoint") +
  ggplot2::scale_x_continuous(breaks = c(1:10), labels = c(1:10)) # this is determined by the order of the FIR

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(file = file.path(graph_fpath_out,
                                "FIR_betas_RHitsVSFhitsMiss_geomline_left-hemi_shaded-errorbars_inclHC.pdf"),
                  width=7, height=6)
}

