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
ROIs <- c("brCA1_body", "brCA2_3_DG_body", "brwhole_hippo")
all_betas <- data.frame()

#' # Load in data
for(isubj in 1:length(subj_formatted)){
  cur_subj <- subj_formatted[isubj]

  for(ihemi in 1:length(hemis)){
    cur_hemi <- hemis[ihemi]

    for(iroi in 1:length(ROIs)){
      cur_roi <- ROIs[iroi]

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
  dplyr::filter(roi %in% c("brCA1_body", "brCA2_3_DG_body")) %>%
  ggplot2::ggplot(ggplot2::aes(cur.mean, fill = roi_lbl)) +
  ggplot2::geom_histogram(breaks = seq(min_val, max_val, by = 2)) +
  ggplot2::facet_grid(subj_id ~ hemi_lbl) +
  ggplot2::xlab("mean beta values") +
  ggplot2::ggtitle("Distribution of beta values") +
  ggplot2::theme(strip.text.y = ggplot2::element_text(size = 10),
                 strip.text.x = ggplot2::element_text(size = 15),
                 axis.text.y = ggplot2::element_text(size = 8))

if(SAVE_GRAPHS_FLAG == 1){
  ggplot2::ggsave(filename = file.path(graph_fpath_out, 'all_subj_mean_beta_hist-overlap_CA1body-CA23DGbody.pdf'),
                  width = 8, height = 10)
}
