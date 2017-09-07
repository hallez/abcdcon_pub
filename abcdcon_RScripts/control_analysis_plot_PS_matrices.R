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
devtools::install_github("rlbarter/superheat")

# add dplyr with library()
# NB: this is non-standard
# correct way would be to make this script into a function,
# store in the R/ directory,
# and use @importFrom dplyr "%>%",
# but this is incompatible with getting in-line results with code in knitr output
library(dplyr)
library(superheat)

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

#' ### Variables needed for plotting loop
subjects <- c("s001", "s002", "s006", "s007", "s008", "s009", "s010",
              "s011", "s012", "s013", "s014", "s018", "s019", "s020",
              "s021", "s023", "s024", "s025", "s026", "s027", "s028", "s029", "s030")
nsub <- length(subjects)
conditions <- c("diffVideo_diffHouse", "sameVideo_sameHouse", "diffVideo_sameHouse")
ncond <- length(conditions)
rois <- c("CA1_body", "CA2_3_DG_body")
nroi <- length(rois)

#' # Define `multiplot` function
# this is from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' # Load in square matrices
# based on `pattern_similarity_no_outlier_trials_load_data_btwn_runs.R`
for(iroi in 1:nroi){

  # reset plotting between ROIs
  allplots <- list()

  for(isubj in 1:nsub){

    cur_subj <- subjects[isubj]
    cur_roi <- rois[iroi]

    cur_file <- file.path(analyzed_mri_dir, cur_subj, "ROIs",
                          "ashs_left", sprintf('br%s_pattern_corr_no_outlier_trials_all_runs.mat',cur_roi))

    if(file.exists(cur_file)){
      # read in the subject's data for the current ROI
      subj_pattern_corr <- data.matrix(as.data.frame(R.matlab::readMat(cur_file)))

      # also read in the labels for the current trials
      subj_pattern_ids <- as.data.frame(R.matlab::readMat(file.path(analyzed_mri_dir, cur_subj, "ROIs",
                                                                    "ashs_left", sprintf('br%s_pattern_mtx_ids_no_outlier_trials_all_runs.mat',cur_roi))))

      # --- format labels so can filter trials ---
      # TODO: figure out why this throws an error (even though output seems reasonable)
      subj_pattern_ids_tidy <- subj_pattern_ids %>%
        tidyr::gather() %>%
        # also split out trial information into separate columns so can easily filter/find indices
        # NB: will get a warning that there are `too many values`, but still separates correctly (I think this is b/c there are multiple fiels for `sep`)
        tidyr::separate(value, into = c("run_trialnum", "trial_lbl"), sep = "Brown|Gray|CR|Miss|FA|Exclude", remove = FALSE) %>%
        dplyr::mutate(runID_trialID = gsub("(^Run)_([01234]*)_(Trial)_([0123456789]*)_", "\\1\\2_\\3\\4", run_trialnum)) %>%
        dplyr::select(-trial_lbl) %>%
        # now pull out trial type information
        tidyr::separate(value, into = c("runlbl", "run_num", "triallbl", "trial_num", "norun_notrial"), sep = "_", remove = FALSE, extra = "merge") %>%
        tidyr::separate(norun_notrial, into = c("norun_notrial_novideo", "video_id"), sep = "EncVideo_", remove = FALSE, extra = "merge") %>%
        tidyr::separate(norun_notrial_novideo, into = c("other_mem_resp", "sep2"), sep = "Brown|Gray", remove = FALSE, extra = "merge") %>%
        tidyr::separate(sep2, into = c("memresp_HI", "memresp_HC"), sep = "Rm1_|Rm2_", remove = FALSE, extra = "merge") %>%
        tidyr::separate(memresp_HC, into = c("memresp_only_HC", "sourceresp_only_HC"), sep = "_", remove = FALSE, extra = "merge") %>%
        # ok, now start combining this information so can easily sort trials
        # item memory
        dplyr::mutate(memresp = ifelse((other_mem_resp!=""), gsub("_", "", other_mem_resp), # NB - this means that false alarms will be `FAF` or `FAR`
                                       ifelse(memresp_HI!="", gsub("_(RHit|FHit)_HI_", "\\1", memresp_HI),
                                              ifelse(!is.na(memresp_only_HC), memresp_only_HC, "999")))) %>%
        dplyr::mutate(memresp_fact = factor(memresp, levels = c("RHit", "FHit", "CR", "FAF", "FAR", "Miss", "ExcludeTrial"))) %>%
        # source memory
        dplyr::mutate(sourceresp = ifelse(video_id == "Lure", NA,
                                          ifelse(memresp_HI!="", gsub("_(RHit|FHit)_(HI)_", "\\2", memresp_HI),
                                                 ifelse(!is.na(sourceresp_only_HC), gsub("_","",sourceresp_only_HC), "999")))) %>%
        # now, create indices for trials that were included in PS
        # remember, this was limited to RHit trials where the house was correctly identified
        dplyr::mutate(mem_PSinclude = ifelse(memresp == "RHit", "yes", "no"),
                      source_PSinclude = ifelse(is.na(sourceresp), "no",
                                                ifelse(sourceresp == "HCRC", "yes",
                                                       ifelse(sourceresp == "HCRI", "yes", "no"))),
                      PS_include = ifelse((mem_PSinclude == "yes" & source_PSinclude == "yes"), "yes", "no"))

      exclude_indices <- which(subj_pattern_ids_tidy$PS_include == "no")

      # save out in case need to use these labels in other analyses
      # technically, labels should be the same for all ROIs w/in a subject,
      # but this is the file naming convention I've used for other analyses so sticking w/ it
      save(subj_pattern_ids_tidy, file = file.path(analyzed_mri_dir, cur_subj, "ROIs",
                                            "ashs_left", sprintf('br%s_pattern_mtx_ids_tidied_no_outlier_trials_all_runs.RData',cur_roi)))

      # --- capitalize on the power of R and give the pattern matrix meaningful row and column names ---
      colnames(subj_pattern_corr) <- subj_pattern_ids_tidy$value
      rownames(subj_pattern_corr) <- subj_pattern_ids_tidy$value

      # --- notch out w/in run autocorrelation ---
      # graphical check: only the NaN-ed out trials should be gray
      # ggplot2::ggplot(reshape2::melt(subj_pattern_corr), ggplot2::aes(Var2, Var1, fill = value)) +
      #   ggplot2::geom_tile() +
      #   ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())
      subj_pattern_corr_no_auto_corr <- subj_pattern_corr
      subj_pattern_corr_no_auto_corr[grep("Run[_]01",rownames(subj_pattern_corr_no_auto_corr)),
                                     grep("Run[_]01",colnames(subj_pattern_corr_no_auto_corr))] <- NA
      subj_pattern_corr_no_auto_corr[grep("Run[_]02",rownames(subj_pattern_corr_no_auto_corr)),
                                     grep("Run[_]02",colnames(subj_pattern_corr_no_auto_corr))] <- NA
      subj_pattern_corr_no_auto_corr[grep("Run[_]03",rownames(subj_pattern_corr_no_auto_corr)),
                                     grep("Run[_]03",colnames(subj_pattern_corr_no_auto_corr))] <- NA
      subj_pattern_corr_no_auto_corr[grep("Run[_]04",rownames(subj_pattern_corr_no_auto_corr)),
                                     grep("Run[_]04",colnames(subj_pattern_corr_no_auto_corr))] <- NA
      # graphical check: now should have gray boxes for each run (where previously could see higher correlations)
      # ggplot2::ggplot(reshape2::melt(subj_pattern_corr_no_auto_corr), ggplot2::aes(Var2, Var1, fill = value)) +
      #   ggplot2::geom_tile() +
      #   ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())

      # --- notch out trials not included in PS analyses ---
      subj_pattern_corr_PS_trials_only <- subj_pattern_corr_no_auto_corr
      subj_pattern_corr_PS_trials_only[exclude_indices, ] <- NA
      subj_pattern_corr_PS_trials_only[,exclude_indices] <- NA
      # graphical check: now should have gray rows and columns for trials not included in PS analyses
      # ggplot2::ggplot(reshape2::melt(subj_pattern_corr_PS_trials_only), ggplot2::aes(Var2, Var1, fill = value)) +
      #   ggplot2::geom_tile() +
      #   ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())

      # --- notch out upper half of symmetric matrix ---
      subj_pattern_corr_lower <- subj_pattern_corr_PS_trials_only
      subj_pattern_corr_lower[lower.tri(subj_pattern_corr_lower)] <- NA
      # graphical check: top of matrix should now all be gray
      # ggplot2::ggplot(reshape2::melt(subj_pattern_corr_lower), ggplot2::aes(Var2, Var1, fill = value)) +
      #   ggplot2::geom_tile() +
      #   ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())

      # --- make pretty plots ---
      if(SAVE_GRAPHS_FLAG == 1){
        png(file.path(graph_fpath_out, sprintf('%s_all-cond_%s_PS_left-hemi_superheat.png', cur_subj, cur_roi)), height = 800, width = 800)
        # can't save plot to a variable b/c it's a list so just graph it while the `png` saver is listening
        superheat::superheat(X = subj_pattern_corr_lower,
                             legend = FALSE,
                             left.label = "none",
                             bottom.label = "none",
                             heat.na.col = "white",
                             heat.lim = c(-1, 1))
        dev.off()
      }

      allplots[[isubj]] <- GGally::ggcorr(subj_pattern_corr_no_auto_corr, size = 0,
                                      legend.position = "none",
                                      low = "#998ec3",
                                      mid = "#f7f7f7",
                                      high = "#f1a340") +
        ggplot2::ggtitle(cur_subj)

      # --- organize plots by video ID ---
      subj_pattern_corr_by_video <- subj_pattern_corr_lower
      # since it's a square matrix, can re-names rows and columns based on IDs
      colnames(subj_pattern_corr_by_video) <- subj_pattern_ids_tidy$video_id
      rownames(subj_pattern_corr_by_video) <- subj_pattern_ids_tidy$video_id
      # ordering based on https://stackoverflow.com/questions/7334644/sort-columns-of-a-dataframe-by-column-name
      subj_pattern_corr_ordered_by_video <- subj_pattern_corr_by_video[, order(colnames(subj_pattern_corr_by_video))]
      subj_pattern_corr_ordered_by_video <- subj_pattern_corr_ordered_by_video[order(rownames(subj_pattern_corr_ordered_by_video)),]
      if(SAVE_GRAPHS_FLAG == 1){
        png(file.path(graph_fpath_out, sprintf('%s_all-cond_%s_PS_left-hemi_by-video_superheat.png', cur_subj, cur_roi)), height = 800, width = 800)
        superheat::superheat(X = subj_pattern_corr_ordered_by_video,
                             legend = FALSE,
                             membership.rows = rownames(subj_pattern_corr_ordered_by_video),
                             left.label.text.size = 3,
                             membership.cols = colnames(subj_pattern_corr_ordered_by_video),
                             bottom.label.text.size = 3,
                             bottom.label.text.angle = 90,
                             heat.na.col = "gray",
                             heat.lim = c(-1, 1))
        dev.off()
      }

      # save out one plot w/ the legend to use for making figures
      if(isubj == 1){
        if(SAVE_GRAPHS_FLAG == 1){
          png(file.path(graph_fpath_out, sprintf('%s_all-cond_%s_PS_left-hemi_by-video_superheat_with-legend.png', cur_subj, cur_roi)), height = 800, width = 800)
          superheat::superheat(X = subj_pattern_corr_ordered_by_video,
                               legend.height = 0.5,
                               legend.width = 2,
                               legend.text.size = 20,
                               membership.rows = rownames(subj_pattern_corr_ordered_by_video),
                               left.label.text.size = 3,
                               membership.cols = colnames(subj_pattern_corr_ordered_by_video),
                               bottom.label.text.size = 3,
                               bottom.label.text.angle = 90,
                               heat.na.col = "gray",
                               heat.lim = c(-1, 1))
          dev.off()
        }
      }

    } #file.exists
  } #isubj
  # save out multiplot before going onto next ROI
  if(SAVE_GRAPHS_FLAG == 1){
    pdf(file.path(graph_fpath_out, sprintf('%s_PS_all-cond_left-hemi_all_subj.pdf', cur_roi)))
    multiplot(plotlist = allplots, cols = 3)
    dev.off()
  }
} #iroi

#' # Load in PS data
# this file is saved out in `mixed_models.R`
load(paste0(analyzed_mri_dir, 'group_z_renamed_spatial_temporal_PS_by_trial.RData'))

#' # Filter out data that's not needed
tidy_trials <- all_trials_z_better_names %>%
  dplyr::filter(hemi == "left",
                condition %in% c("diffVideo_diffHouse", "sameVideo_sameHouse", "diffVideo_sameHouse"),
                roi %in% c("CA1_body", "CA2_3_DG_body")) %>%
  dplyr::select(-z_r)

#' ## Compute subject-wise means for each condition and ROI
mean_vals <- tidy_trials %>%
  dplyr::group_by(subj, roi, condition) %>%
  dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
  tidyr::spread(condition, mean_r) %>%
  # revalue conditions so fit w/ other plots in the paper
  dplyr::rename("Different Video\nSame House" = diffVideo_sameHouse,
                "Different Video\nDifferent House" = diffVideo_diffHouse,
                "Same Video\nSame House" = sameVideo_sameHouse)

#' ## Plot by ROI
for(iroi in 1:nroi){
  cur_roi <- rois[iroi]
  cur_mean_vals <- mean_vals %>%
    dplyr::ungroup() %>%
    dplyr::filter(roi == cur_roi) %>%
    as.data.frame()

  # `dplyr::add_rownames` and `tibble::column_to_rownames` are both depricated, so do this instead
  mean_vals_rowids <- cur_mean_vals
  rownames(mean_vals_rowids) <- mean_vals_rowids$subj
  mean_vals_by_subj <- mean_vals_rowids %>%
    dplyr::select(-subj, -roi)

  if(SAVE_GRAPHS_FLAG == 1){
    png(file.path(graph_fpath_out, sprintf('all-cond-means_%s_PS_left-hemi.png', cur_roi)), height = 800, width = 800)
    superheat::superheat(X = mean_vals_by_subj,
                         title = gsub("_body","",cur_roi))
    dev.off()
  }

}

#' ## plot both ROIs in same plot so that scales match
mean_vals_by_ROI <- tidy_trials %>%
  dplyr::group_by(subj, roi, condition) %>%
  dplyr::summarise(mean_r = mean(r, na.rm = TRUE)) %>%
  tidyr::unite(condition_roi, condition, roi) %>%
  tidyr::spread(condition_roi, mean_r) %>%
  # revalue conditions so fit w/ other plots in the paper
  dplyr::rename("CA1\nDifferent Video\nSame House" = diffVideo_sameHouse_CA1_body,
                "CA1\nDifferent Video\nDifferent House" = diffVideo_diffHouse_CA1_body,
                "CA1\nSame Video\nSame House" = sameVideo_sameHouse_CA1_body,
                "CA23DG\nDifferent Video\nSame House" = diffVideo_sameHouse_CA2_3_DG_body,
                "CA23DG\nDifferent Video\nDifferent House" = diffVideo_diffHouse_CA2_3_DG_body,
                "CA23DG\nSame Video\nSame House" = sameVideo_sameHouse_CA2_3_DG_body) %>%
  dplyr::ungroup() %>%
  as.data.frame()

mean_vals_by_ROI_rowids <- mean_vals_by_ROI
rownames(mean_vals_by_ROI_rowids) <- mean_vals_by_ROI_rowids$subj
mean_vals_by_ROI_by_subj <- mean_vals_by_ROI_rowids %>%
  dplyr::select(-subj)

if(SAVE_GRAPHS_FLAG == 1){
  png(file.path(graph_fpath_out, sprintf('all-cond-means_CA1-CA23DG_columns-by-cond_PS_left-hemi.png')), height = 800, width = 800)
  superheat::superheat(X = mean_vals_by_ROI_by_subj)
  dev.off()
}

# reorder by ROI
mean_vals_by_ROI_by_subj_order_by_ROI <- mean_vals_by_ROI_by_subj %>%
  dplyr::select(starts_with("CA1"), starts_with("CA23DG"))

if(SAVE_GRAPHS_FLAG == 1){
  png(file.path(graph_fpath_out, sprintf('all-cond-means_CA1-CA23DG_columns-by-ROI_PS_left-hemi.png')), height = 800, width = 800)
  superheat::superheat(X = mean_vals_by_ROI_by_subj_order_by_ROI)
  dev.off()
}

#' ## Loop through and plot from PS data
for(iroi in 1:nroi){
  for(icond in 1:ncond){
    # right now, save out all plots for a given condition together
    # based on: https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
    myplots <- list()

    for(isubj in 1:nsub){
      cur_subj <- subjects[isubj]
      cur_cond <- conditions[icond]
      cur_roi <- rois[iroi]

      # --- plot just the current condition's trial pairs ---
      cur_dat <- tidy_trials %>%
        dplyr::filter(subj == cur_subj,
                      condition == cur_cond,
                      roi == cur_roi)

      cur_dat_fmt <- cur_dat %>%
        tidyr::spread(col_name, r) %>%
        dplyr::select(-subj, -roi, -hemi, -condition) %>%
        # based on: https://stackoverflow.com/questions/5555408/convert-the-values-in-a-column-into-row-names-in-an-existing-data-frame-in-r
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = "row_name")

      myplots[[isubj]] <- GGally::ggcorr(cur_dat_fmt, size = 0,
                                         legend.position = "none",
                                         low = "#998ec3",
                                         mid = "#f7f7f7",
                                         high = "#f1a340")

      if(SAVE_GRAPHS_FLAG == 1){
        ggplot2::ggsave(file = file.path(graph_fpath_out, sprintf('%s_%s_%s_PS_left-hemi.pdf', cur_subj, cur_cond, cur_roi)),
                        height = 5, width = 5)
      }

      # also plot a heatmap version
      # if want to match colors to `ggcorr` - heat.pal = c("#998ec3", "white", "#f1a340")
      lower_dat <- cur_dat_fmt
      lower_dat[lower.tri(lower_dat)] <- NA
      if(SAVE_GRAPHS_FLAG == 1){
        png(file.path(graph_fpath_out, sprintf('%s_%s_%s_PS_left-hemi_superheat.png', cur_subj, cur_cond, cur_roi)), height = 800, width = 800)
        # can't save plot to a variable b/c it's a list so just graph it while the `png` saver is listening
        superheat::superheat(X = lower_dat,
                             legend = FALSE,
                             left.label = "none",
                             bottom.label = "none",
                             heat.na.col = "white")
        dev.off()
      }
    } #isubj

    # in order to save a multiplot, need to have the file device that will be saving to open
    # based on: https://stackoverflow.com/questions/11721401/r-save-multiplot-to-file
    if(SAVE_GRAPHS_FLAG == 1){
      pdf(file.path(graph_fpath_out, sprintf('%s_%s_PS_left-hemi_all_subj.pdf', cur_cond, cur_roi)))
      multiplot(plotlist = myplots, cols = 2)
      dev.off()
    }

  } #icond
} #iroi

#' # Use `ggplot2::geom_tile` for plotting
# based on: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
tidy_trials %>%
  # this seems like it doesn't work, but it's just struggling w/ the scaling across so many subjects
  # uncomment this out to see it work:
  # dplyr::filter(subj %in% c("s001", "s002", "s018") %>%
  dplyr::filter(hemi == "left",
                roi == "CA1_body",
                condition == "diffVideo_diffHouse") %>%
  ggplot2::ggplot(ggplot2::aes(x = row_name, y = col_name, fill = r)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank()) +
  ggplot2::facet_grid(subj ~ ., scales = "free")
