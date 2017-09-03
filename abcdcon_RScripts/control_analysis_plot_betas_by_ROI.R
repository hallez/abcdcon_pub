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
                          halle::ensure_trailing_slash("figures"),
                          halle::ensure_trailing_slash("plot-beta-distributions"))

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
all_betas <- data.frame()

#' # Load in data
for(isubj in 1:length(subj_formatted)){
  cur_subj <- subj_formatted[isubj]

  for(ihemi in 1:length(hemis)){
    cur_hemi <- hemis[ihemi]

    for(iroi in 1:length(ROIs)){
      cur_roi <- ROIs[iroi]

      # this file is created by `control_analysis_betas_by_ROI.m`
      cur_betas_fpath <- file.path(graph_fpath_out, sprintf('%s_%s_%s_means.mat', cur_subj, cur_hemi, cur_roi))

      if(file.exists(cur_betas_fpath)){
        cur_betas <- as.data.frame(R.matlab::readMat(cur_betas_fpath))

        # --- Add some labels to the `cur_betas` dataframe ---
        cur_betas_tidy <- cur_betas %>%
          dplyr::mutate(subj_id = cur_subj,
                        hemi = cur_hemi,
                        roi = cur_roi,
                        voxel_id = row_number())

        # --- Merge into group dataframe ---
        if(dim(all_betas)[1]==0){
          all_betas <- cur_betas_tidy
        } else {
          all_betas <- dplyr::full_join(all_betas, cur_betas_tidy, by = intersect(names(all_betas), names(cur_betas_tidy)))
        } #dim(all_betas
      } else {
        print(sprintf("Current mean betas file %s does not exist. Continuing on.", cur_betas_fpath))
        next
      } #if(file.exists
    } #iroi
  } #ihemi
} #isubj

#' # Tidy up group dataframe
#' ## Print out untiedied info
head(all_betas)
unique(all_betas$subj_id)

#' ## Tidy up labels
all_betas_tidy <-
  all_betas %>%
  dplyr::mutate(hemi_lbl = gsub("ashs_", "", hemi),
                roi_lbl = gsub("br","",roi))
#' ## Peek at tidied data
head(all_betas_tidy)

#' # Test whether mean beta value differs by ROI
#' ## Mean
mean_betas <- all_betas_tidy %>%
  dplyr::group_by(subj_id, roi) %>%
  dplyr::summarise(mean_val = mean(cur.mean, na.rm = TRUE)) %>%
  tidyr::spread(key = roi, value = mean_val)

t.test(mean_betas$brCA1_body, mean_betas$brCA2_3_DG_body, paired = TRUE)

mean_betas %>%
  dplyr::ungroup() %>%
  dplyr::rename(var1 = brCA1_body, var2 = brCA2_3_DG_body) %>%
  dplyr::select(var1, var2) %>%
  halle::compute_cohens_d()

#' ## Median
median_betas <- all_betas_tidy %>%
  dplyr::group_by(subj_id, roi) %>%
  dplyr::summarise(median_val = median(cur.mean, na.rm = TRUE)) %>%
  tidyr::spread(key = roi, value = median_val)

t.test(median_betas$brCA1_body, median_betas$brCA2_3_DG_body, paired = TRUE)

median_betas %>%
  dplyr::ungroup() %>%
  dplyr::rename(var1 = brCA1_body, var2 = brCA2_3_DG_body) %>%
  dplyr::select(var1, var2) %>%
  halle::compute_cohens_d()

#' # Plot
#' ## Figure out range of betas to set limits
max_val <- ceiling(max(all_betas_tidy$cur.mean, na.rm = TRUE))
min_val <- floor(min(all_betas_tidy$cur.mean, na.rm = TRUE))

#' ## Define plotting functions
plot_hist <- function (df_in, roi_val) {
  p <- NULL

  p <<- df_in %>%
    ggplot2::ggplot(ggplot2::aes(cur.mean)) +
    ggplot2::geom_histogram(breaks = seq(min_val, max_val, by = 2)) +
    ggplot2::facet_grid(subj_id ~ hemi_lbl) +
    ggplot2::xlab("mean beta values") +
    ggplot2::ggtitle(sprintf("Distribution of beta values for %s", roi_val)) +
    ggplot2::theme(strip.text.y = ggplot2::element_text(size = 10),
                   strip.text.x = ggplot2::element_text(size = 15),
                   axis.text.y = ggplot2::element_text(size = 8))
}

#' ## Plot by ROI
for(iroi in 1:length(ROIs)){
  cur_roi <- ROIs[iroi]
  roi_lbl <- unique(all_betas_tidy$roi_lbl[all_betas_tidy$roi == cur_roi])

  all_betas_tidy %>%
    dplyr::filter(roi == cur_roi) %>%
    plot_hist(., roi_lbl)

  print(p)

  if(SAVE_GRAPHS_FLAG == 1){
    ggplot2::ggsave(filename = file.path(graph_fpath_out, sprintf('all_subj_mean_beta_hist_%s.pdf', cur_roi)),
                    width = 8, height = 10)
  }
}

#' ## Overlay distibutions for CA1 and CA23DG (body)
all_betas_tidy %>%
  dplyr::filter(roi %in% c("brCA1_body", "brCA2_3_DG_body"),
                hemi == "ashs_left") %>%
  # format ROI labels to match other plots
  dplyr::mutate(roi_lbl = gsub("(.*?)_body", "\\1", roi_lbl),
                roi_lbl = sub("_","",roi_lbl),
                roi_lbl = sub("_","",roi_lbl)) %>%
  ggplot2::ggplot(ggplot2::aes(color = roi_lbl, fill = roi_lbl)) +
  ggplot2::geom_violin(ggplot2::aes(x = roi_lbl, y = cur.mean), alpha = 0.6) +
  ggplot2::facet_wrap(~subj_id) +
  ggplot2::ylab("Beta Value") +
  ggplot2::theme(strip.text.y = ggplot2::element_text(size = 10),
                 axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8),
                 legend.title = ggplot2::element_blank())

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(filename = file.path(graph_fpath_out, 'all_subj_mean_beta_violin_CA1body-CA23DGbody.pdf'),
                  width = 8, height = 10)
}

#' # Plot normalizing for the number of voxels
# w/ a single ROI, this will scale to 1
all_betas_tidy %>%
  dplyr::filter(roi == "brCA1_body") %>%
  ggplot2::ggplot(ggplot2::aes(cur.mean, fill = roi_lbl)) +
  ggplot2::geom_histogram(ggplot2::aes(y = ..ncount..), breaks = seq(min_val, max_val, by = 2)) +
  ggplot2::facet_grid(subj_id ~ hemi_lbl) +
  ggplot2::xlab("mean beta values") +
  ggplot2::theme(strip.text.y = ggplot2::element_text(size = 10),
                 strip.text.x = ggplot2::element_text(size = 15),
                 axis.text.y = ggplot2::element_text(size = 8))

# w/ 2 ROIs, this scales to 2
all_betas_tidy %>%
  dplyr::filter(roi %in% c("brCA1_body", "brCA2_3_DG_body")) %>%
  ggplot2::ggplot(ggplot2::aes(cur.mean, fill = roi_lbl)) +
  ggplot2::geom_histogram(ggplot2::aes(y = ..ncount..), breaks = seq(min_val, max_val, by = 2)) +
  ggplot2::facet_grid(subj_id ~ hemi_lbl) +
  ggplot2::xlab("mean beta values") +
  ggplot2::theme(strip.text.y = ggplot2::element_text(size = 10),
                 strip.text.x = ggplot2::element_text(size = 15),
                 axis.text.y = ggplot2::element_text(size = 8))

# need to use `..density..` to get it to scale for each group (but then y-axis isn't out of 1)
all_betas_tidy %>%
  dplyr::filter(roi %in% c("brCA1_body", "brCA2_3_DG_body")) %>%
  # normalized yaxis based on: https://stackoverflow.com/questions/11766856/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion
  # use `..density..` b/c then it's normalized by group (see: https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group)
  # need to scale by `binwidth` b/c this will scale it appropriately
  ggplot2::ggplot(ggplot2::aes(x = cur.mean, fill = roi_lbl)) +
  ggplot2::geom_histogram(ggplot2::aes(y = 2*..density..), breaks = seq(min_val, max_val, by = 2), binwidth = 2) +
  ggplot2::facet_grid(subj_id ~ hemi_lbl) +
  ggplot2::xlab("mean beta values") +
  ggplot2::ylab("density") +
  ggplot2::ggtitle("Distribution of beta values") +
  ggplot2::theme(strip.text.y = ggplot2::element_text(size = 10),
                 strip.text.x = ggplot2::element_text(size = 15),
                 axis.text.y = ggplot2::element_text(size = 8))
