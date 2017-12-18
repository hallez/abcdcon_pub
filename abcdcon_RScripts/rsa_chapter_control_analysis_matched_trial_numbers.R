#' ---
#' title: ABCDCon Matched Trial Numbers
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
graph_fpath_out <- file.path("Users", "hrzucker", "Dropbox", "work", "manuscripts", "rsa-methods-chapter", "figures")

#' ## Setup other variables
#' ### Flags
SAVE_GRAPHS_FLAG <- 1

#+ label="Load data"
num_reps <- 1000
load(file = paste0(analyzed_mri_dir, sprintf('match_trialnums_%dreps_chisqvals.RData', num_reps)))

#' # Plot histogram of chi-square values
all_chisq %>%
  ggplot2::ggplot(ggplot2::aes(chisq_value)) +
  ggplot2::geom_histogram(binwidth = 0.25)

ggplot2::ggsave(file = paste0(graph_fpath_out,
                              halle::ensure_trailing_slash("matched-trial-nums"),
                              sprintf("match_trialnums_%dreps_chisqvals_histogram.pdf", num_reps)),
                width=8, height=6)

#' # Graph
#+ echo = FALSE
group_subset %>%
  dplyr::filter(condition != "anyVideo_sameHouse") %>%
  dplyr::filter(hemi == "left") %>%
  dplyr::group_by(roi, hemi, condition) %>%
  dplyr::summarise(mean = mean(r),
                   sd = sd(r),
                   n = length(r),
                   sem = sd(r)/sqrt(length(r))) %>%
  # re-order conditions
  dplyr::mutate(condition_ordered = factor(condition, levels = c("sameVideo_sameHouse", "diffVideo_sameHouse", "diffVideo_diffHouse"))) %>%
  # put one space between same/diff and video/house and TWO spaces between where we'll want a line break
  # this will help `gsub` only create 2 lines for each condition name rather than 4 lines
  dplyr::mutate(condition_renamed = car::recode(condition_ordered, "'diffVideo_sameHouse' = 'Different Video  Same House'; 'diffVideo_diffHouse' = 'Different Video  Different House' ; 'sameVideo_sameHouse' = 'Same Video  Same House'")) %>%
  dplyr::mutate(condition_renamed_ordered = factor(condition_renamed, levels = c("Same Video  Same House", "Different Video  Same House", "Different Video  Different House"))) %>%
  # based on https://www.r-bloggers.com/line-breaks-between-words-in-axis-labels-in-ggplot-in-r/
  dplyr::mutate(condition_breaks = gsub("  ", "\n", condition_renamed_ordered)) %>%
  ggplot2::ggplot(ggplot2::aes(x = condition_breaks, y = mean, fill = condition_breaks)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::facet_grid(. ~ roi) +
  # to match colors from method figure:
  # "#CC6633" = orange (diff video, diff house)
  # "#CC6699" = fuscia (same video, same house)
  # "#66CC33" = green (different video, same house)
  ggplot2::scale_fill_manual(values = c("#CC6633", "#66CC33", "#CC6699")) +
  ggplot2::geom_errorbar(ggplot2::aes(ymax = mean + sem,
                                      ymin = mean - sem,
                                      width=0.10)) +
  ggplot2::ggtitle("Neural pattern similarity for spatial and event contexts") +
  ggplot2::ylab("Mean Pattern Similarity (r)") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, color = "black"), axis.title.x = ggplot2::element_blank(),
                 strip.text.x = ggplot2::element_text(size = 20),
                 axis.text.y = ggplot2::element_text(size = 10), axis.title.y = ggplot2::element_text(size = 20),
                 strip.text.y = ggplot2::element_text(size = 20),
                 legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=20, vjust=2)) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(file = paste0(graph_fpath_out,
                              halle::ensure_trailing_slash("matched-trial-nums"),
                              "sameAll_sameSome_diffAll_left_body_CA1_CA23DG_matched_trial_nums.pdf"),
                width=8, height=6)


