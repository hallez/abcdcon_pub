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
dropbox_graph_fpath_out <- paste0(halle::ensure_trailing_slash(dropbox_dir),
                          halle::ensure_trailing_slash("writeups"),
                          halle::ensure_trailing_slash("figures"))
graph_fpath_out <- paste0(dropbox_graph_fpath_out, "RT-PS_correlation")
dir.create(graph_fpath_out) # will throw an error if this already exists

#' ## Setup other variables
rois_of_interest <- c("CA1_body", "CA2_3_DG_body")
GRAPH_FLAG <- 1
SAVE_GRAPHS_FLAG <-1
MODEL_FLAG <- 1

#' # Load in PS data
load(file.path(analyzed_mri_dir, 'group_z_renamed_spatial_temporal_PS_by_trial.RData'))

#' ## Peek at the PS data
head(all_trials_z_better_names)
unique(all_trials_z_better_names$roi)
unique(all_trials_z_better_names$condition)
unique(all_trials_z_better_names$subj)

#' ## Setup conditions of interest
conditions_of_no_interest <- "anyVideo_sameHouse"
conditions_of_interest <- unique(all_trials_z_better_names$condition)
conditions_of_interest <- conditions_of_interest[conditions_of_interest != conditions_of_no_interest]

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

#' ## Peek at merged data
head(alldat)

#' ## Sanity spot check RTs
#' ### Rows
alldat %>%
  dplyr::filter(row.runID_trialID == "Run02_Trial010",
                hemi == "left",
                roi == "CA1_body",
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
                roi == "CA1_body",
                subj == "s006") %>%
  dplyr::distinct(MemRT_col)

behav_trim %>%
  dplyr::filter(col.runID_trialID == "Run04_Trial038",
                subj == "s006") %>%
  dplyr::distinct(MemRT)

#' # Compute mean (and SD) RT difference by PS condition (Supplemental Table 1)
alldat %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(mean_RT_diff = mean(RT_diff, na.rm = TRUE),
                   sd_RT_diff = sd(RT_diff, na.rm = TRUE))

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
# limit this to left hemi since this is where we find sigificant effects
if(GRAPH_FLAG == 1) {
  allplots <- list()
  plot_counter <- 0
  for(iroi in 1:length(rois_of_interest)){
    for(icond in 1:length(conditions_of_interest)){
      cur_roi <- rois_of_interest[iroi]
      cur_cond <- conditions_of_interest[icond]

      # --- plot all subjects together ---
      alldat %>%
        dplyr::filter(hemi == "left",
                      roi == cur_roi,
                      condition == cur_cond) %>%
        scatter_cond()
      p <- p +
        ggplot2::ggtitle(sprintf("%s: %s", cur_roi, cur_cond)) +
        ggplot2::xlab("trial pair correlation value (r)") +
        ggplot2::ylab("trial pair RT difference (in ms)")
      print(p)

      if(SAVE_GRAPHS_FLAG == 1){
        ggplot2::ggsave(file = file.path(graph_fpath_out,
                                         sprintf("RT-PS_correlation_%s_%s.pdf", cur_roi, cur_cond)),
                        width=8, height=6)
      }

      # --- plot faceting by subject ---
      alldat %>%
        dplyr::filter(hemi == "left",
                      roi == cur_roi,
                      condition == cur_cond) %>%
        scatter_cond_facet_by_subj()
      p <- p +
        ggplot2::ggtitle(sprintf("%s: %s \n by subject ", cur_roi, cur_cond)) +
        ggplot2::xlab("trial pair correlation value (r)") +
        ggplot2::ylab("trial pair RT difference (in ms)") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6))
      print(p)


      if(SAVE_GRAPHS_FLAG == 1){
        ggplot2::ggsave(file = file.path(graph_fpath_out,
                                         sprintf("RT-PS_correlation_%s_%s_bySubj.pdf", cur_roi, cur_cond)),
                        width=8, height=6)

        # --- format for figure ---
        plot_counter <- plot_counter + 1
        cur_roi_fmt <- gsub("(.*?)_body", "\\1",cur_roi) %>%
          sub("_","",.) %>%
          sub("_","",.)
        cur_cond_fmt <- ifelse((cur_cond == "diffVideo_sameHouse"), "Different Video Same House",
                               ifelse((cur_cond == "diffVideo_diffHouse"), "Different Video Different House",
                               ifelse((cur_cond == "sameVideo_sameHouse"), "Same Video Same House", "error")))

        cp <- p +
          # "pretty" scale from: https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
          ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
          ggplot2::ggtitle(sprintf("%s: %s", cur_roi_fmt, cur_cond_fmt)) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6),
                         # skinny margins from: https://stackoverflow.com/questions/36783189/changing-the-appearance-of-facet-labels-size
                         strip.text.x = ggplot2::element_text(size = 10, margin = margin(.1, 0, .1, 0, "cm")))
        allplots[[plot_counter]] <- cp

      }
    } #icond
  } #iroi
  plot3x2 <- cowplot::plot_grid(allplots[[1]], allplots[[2]], allplots[[3]], allplots[[4]], allplots[[5]], allplots[[6]],
                     labels = c("A.", "B.", "C.", "D.", "E.", "F."),
                     ncol = 3, nrow = 2)
  cowplot::save_plot(file.path(graph_fpath_out, "RT-PS_correlation_all-cond_all-roi_bySubj.png"),
            plot3x2,
            ncol = 3,
            nrow = 2,
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3
  )

} # GRAPH_FLAG

#' # Re-model including RT
if(MODEL_FLAG == 1) {
  #' ## Setup dataframe
  all_trials_body <-
    alldat %>%
    dplyr::filter(roi %in% c("CA1_body", "CA2_3_DG_body")) %>%
    droplevels(.)

  head(all_trials_body)
  unique(all_trials_body$roi)

  #' ## filter so just have conditions of interest
  all_trials_body_SVSH_DVSH_DVDH <- NULL
  all_trials_body_SVSH_DVSH_DVDH <-
    all_trials_body %>%
    dplyr::filter(condition %in% c("sameVideo_sameHouse", "diffVideo_sameHouse", "diffVideo_diffHouse"))
  head(all_trials_body_SVSH_DVSH_DVDH)
  unique(all_trials_body_SVSH_DVSH_DVDH$condition)
  unique(all_trials_body_SVSH_DVSH_DVDH$hemi)
  unique(all_trials_body_SVSH_DVSH_DVDH$roi)

  #' ### main effects w/o random slopes
  lmm.all_cond_roi_hemi.no_random_slopes <- lme4::lmer(z_r ~ condition*roi + condition*hemi + roi*hemi + (1|RT_diff) + (1|subj), data = all_trials_body_SVSH_DVSH_DVDH, REML = FALSE)
  summary(lmm.all_cond_roi_hemi.no_random_slopes)

  #' ### 3-way interaction w/o random slopes
  lmm.all_condXroiXhemi.no_random_slopes <- lme4::lmer(z_r ~ condition*roi*hemi + (1|RT_diff) + (1|subj), data = all_trials_body_SVSH_DVSH_DVDH, REML = FALSE)
  summary(lmm.all_condXroiXhemi.no_random_slopes)

  #' #### Model comparisons (ROI x Context Similarity x Hemisphere)
  print(cat("\nModel comparisons (ROI x Context Similarity x Hemisphere) including (1|RT_diff)"))
  print(anova(lmm.all_cond_roi_hemi.no_random_slopes, lmm.all_condXroiXhemi.no_random_slopes))

  #' # Left hemi only
  all_trials_body_SVSH_DVSH_DVDH_left <-
    all_trials_body_SVSH_DVSH_DVDH %>%
    dplyr::filter(hemi == "left")

  lmm.all_cond_roi.left <- lme4::lmer(z_r ~ condition + roi + (1|RT_diff) + (1|subj), data = all_trials_body_SVSH_DVSH_DVDH_left, REML = FALSE)
  summary(lmm.all_cond_roi.left)

  lmm.all_condXroi.left <- lme4::lmer(z_r ~ condition*roi + (1|RT_diff) + (1|subj), data = all_trials_body_SVSH_DVSH_DVDH_left, REML=FALSE)
  summary(lmm.all_condXroi.left)

  #' ## model comparisons (Left Hemi: ROI x Context Similarity)
  print(cat("\nModel comparisons (Left Hemi: ROI x Context Similarity) including (1|RT_diff)"))
  print(anova(lmm.all_cond_roi.left, lmm.all_condXroi.left))

  #' # Test for interactions between conditions - Episodic Context Similarity
  #' ## Print out what's in dataframe before start running stats
  # yes, this is redundant with when the dataframe gets setup
  # but better to be redundant and ensure you know what you're working with!
  all_trials_body_temporal <- NULL
  all_trials_body_temporal <-
    all_trials_body %>%
    dplyr::filter(condition %in% c("sameVideo_sameHouse", "diffVideo_sameHouse"))
  head(all_trials_body_temporal)
  unique(all_trials_body_temporal$condition)
  unique(all_trials_body_temporal$hemi)
  unique(all_trials_body_temporal$roi)

  #' ## Model setup: Left
  # Now, we need to follow up the significant 3-way interaction in left hemi
  temporal_lmm_left <-
    all_trials_body_temporal %>%
    dplyr::filter(hemi == "left")

  # print out what's in the dataframe so we're supersure before running stats
  head(temporal_lmm_left)
  unique(temporal_lmm_left$condition)
  unique(temporal_lmm_left$roi)
  unique(temporal_lmm_left$hemi)

  # condition, roi
  lmm.ss2 <- lme4::lmer(z_r ~ condition + roi + (1|RT_diff) + (1|subj), data = temporal_lmm_left, REML = FALSE)
  summary(lmm.ss2)

  # condition * roi
  lmm.ss3 <- lme4::lmer(z_r ~ condition*roi + (1|RT_diff) + (1|subj), data = temporal_lmm_left, REML=FALSE)
  summary(lmm.ss3)

  #' ### Compare models (left CA1 x left CA23DG: Episodic Context Similarity)
  # roi and condition vs. roi*condition
  # this is the most comparable analysis to the condition x roi x hemi interaction model we're trying to breakdown
  print(cat("\nCompare models (left CA1 x left CA23DG: Episodic Context Similarity) including (1|RT_diff)"))
  print(anova(lmm.ss2, lmm.ss3))

  #' # Test for interactions between conditions - Spatial Context Similarity
  #' ## Print out what's in dataframe before start running stats
  # yes, this is redundant with when the dataframe gets setup
  # but better to be redundant and ensure you know what you're working with!
  all_trials_body_spatial <- NULL
  all_trials_body_spatial <-
    all_trials_body %>%
    dplyr::filter(roi %in% c("CA1_body", "CA2_3_DG_body")) %>%
    dplyr::filter(condition %in% c("diffVideo_sameHouse", "diffVideo_diffHouse"))
  head(all_trials_body_spatial)
  unique(all_trials_body_spatial$condition)
  unique(all_trials_body_spatial$hemi)
  unique(all_trials_body_spatial$roi)

  #' ## Model setup: Left
  # Now, we need to follow up the significant 3-way interaction in left hemi
  spatial_lmm_left <-
    all_trials_body_spatial %>%
    dplyr::filter(hemi == "left")

  # print out what's in the dataframe so we're supersure before running stats
  head(spatial_lmm_left)
  unique(spatial_lmm_left$condition)
  unique(spatial_lmm_left$roi)
  unique(spatial_lmm_left$hemi)

  # condition, roi
  lmm.ss2 <- lme4::lmer(z_r ~ condition + roi + (1|RT_diff) + (1|subj), data = spatial_lmm_left, REML = FALSE)
  summary(lmm.ss2)

  # condition * roi
  lmm.ss3 <- lme4::lmer(z_r ~ condition*roi + (1|RT_diff) +(1|subj), data = spatial_lmm_left, REML=FALSE)
  summary(lmm.ss3)

  #' ### Compare models (left CA1 x left CA23DG: Spatial Context Similarity)
  # roi and condition vs. roi*condition
  # this is the most comparable analysis to the condition x roi x hemi interaction model we're trying to breakdown
  print(cat("\nCompare models (left CA1 x left CA23DG: Spatial Context Similarity) (1|RT_diff)"))
  print(anova(lmm.ss2, lmm.ss3))

} #MODEL_FLAG
