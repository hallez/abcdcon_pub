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

#+ label="Load data"
#' ## Load in PS data
# this file is saved out in `mixed_models.R`
load(paste0(analyzed_mri_dir, 'group_z_renamed_spatial_temporal_PS_by_trial.RData'))

#' # Filter out data that's not needed
tidy_trials <- all_trials_z_better_names %>%
  dplyr::filter(hemi == "left",
                condition %in% c("diffVideo_diffHouse", "sameVideo_sameHouse", "diffVideo_sameHouse"),
                roi %in% c("CA1_body", "CA2_3_DG_body")) %>%
  dplyr::select(-z_r)

#' # Setup variables needed for plotting loop
subjects <- unique(tidy_trials$subj)
nsub <- length(subjects)
conditions <- unique(tidy_trials$condition)
ncond <- length(conditions)
rois <- unique(tidy_trials$roi)
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

#' # Loop through and plot
for(iroi in 1:nroi){
  for(icond in 1:ncond){
    # right now, save out all plots for a given condition together
    # based on: https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
    myplots <- list()

    for(isubj in 1:nsub){
      cur_subj <- subjects[isubj]
      cur_cond <- conditions[icond]
      cur_roi <- rois[iroi]

      cur_dat <- tidy_trials %>%
        dplyr::filter(subj == cur_subj,
                      condition == cur_cond,
                      roi == cur_roi)

      cur_dat_fmt <- cur_dat %>%
        tidyr::spread(col_name, r) %>%
        dplyr::select(-subj, -roi, -hemi, -condition, -row_name)

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
      p <- superheat::superheat(lower_dat,
                           legend = FALSE,
                           left.label = "none",
                           bottom.label = "none",
                           heat.na.col = "white")
      if(SAVE_GRAPHS_FLAG == 1){
        png(file.path(graph_fpath_out, sprintf('%s_%s_%s_PS_left-hemi_superheat.png', cur_subj, cur_cond, cur_roi)), height = 900, width = 800)
        superheat(X = lower_dat, scale = F)
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
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank()) +
  ggplot2::facet_grid(subj ~ ., scales = "free")
