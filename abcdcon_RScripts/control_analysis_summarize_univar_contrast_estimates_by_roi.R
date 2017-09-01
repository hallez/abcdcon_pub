#' ---
#' title: ABCDCon Univariate Within ROI Estimates
#' author: Halle R. Dimsdale-Zucker
#' output:
#'  html_document:
#'    toc: true
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
SAVE_GRAPHS_FLAG <- 1

#' ## Set paths as variables
analyzed_mri_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_mri))
raw_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$raw_behavioral))
analyzed_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_behavioral))
contrast_est_path <- paste0(analyzed_mri_dir,halle::ensure_trailing_slash('univariate_sanityCheck'))
contrast_est_fname <- 'cons_body_ROIs.csv' # this is created by `control_analysis_univariate_contrast_estimates_by_roi.m`
contrast_est_fpath <- paste0(contrast_est_path,contrast_est_fname)
dropbox_dir <- halle::ensure_trailing_slash(config$directories$dropbox_abcdcon)
dropbox_graph_fpath_out <- paste0(halle::ensure_trailing_slash(dropbox_dir),
                                  halle::ensure_trailing_slash("writeups"),
                                  halle::ensure_trailing_slash("figures"))
graph_fpath_out <- paste0(dropbox_graph_fpath_out, "univariate-by-subfield")
dir.create(graph_fpath_out) # will throw an error if this already exists

if(file.exists(contrast_est_fpath)){
  cons <- NULL
  cons <- read.csv(contrast_est_fpath)

  signif_level <- 0.05/5 # this is bonferroni correction (one hemi)

  cons %>%
    # Hadley-ify formatting
    # one row per observation
    tidyr::spread(contrast,activity) %>%
    # remove subject column because it's meaningless
    dplyr::select(-subject) %>%
    dplyr::group_by(roi_file,hemi) %>%
    dplyr::summarise_each(funs(t.test(.,mu=0,na.rm=TRUE)$p.value)) %>%
    dplyr::group_by(roi_file,hemi) %>%
    dplyr::summarise_each(funs(signif(.,3))) %>%
    dplyr::group_by(roi_file,hemi) %>%
    # get snazzy and add * if signif
    dplyr::summarise_each(funs(ifelse(.<signif_level,paste0(sprintf("%.3f", .),"*"),sprintf("%.3f", .))))

  # --- Rhits vs Fhits & misses ---
  cons %>%
    # Hadley-ify formatting
    # one row per observation
    tidyr::spread(contrast,activity) %>%
    # remove subject column because it's meaningless
    tidyr::gather(condition, value, -roi_file,-hemi,-subject) %>%
    dplyr::group_by(roi_file,hemi,condition) %>%
    dplyr::summarize(activity=mean(value),
              sem=sd(value)/sqrt(length(value)),
              lower = activity - sem,
              upper = activity + sem) %>%
    dplyr::filter(roi_file=="rCA1_body" | roi_file=="rCA2_3_DG_body" | roi_file=="rwhole_hippo",
           condition == "allRHits_vs_FamHitsANDMiss") %>%
    ggplot2::ggplot(ggplot2::aes(condition,activity,fill=roi_file)) +
    ggplot2::geom_bar(width=0.7,position=ggplot2::position_dodge(0.9),stat="identity") +
    ggplot2::facet_grid(.~hemi) +
    ggplot2::ylab("mean activity") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), strip.text.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 20), axis.title.y=ggplot2::element_text(size = 20),
          legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size=20)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lower,ymax=upper),position=ggplot2::position_dodge(width=0.9),color="black",width=0.25)

  # --- split up hits by house ---
  cons %>%
    # Hadley-ify formatting
    # one row per observation
    tidyr::spread(contrast,activity) %>%
    # remove subject column because it's meaningless
    tidyr::gather(condition, value, -roi_file,-hemi,-subject) %>%
    dplyr::group_by(roi_file,hemi,condition) %>%
    dplyr::summarize(activity=mean(value),
                     sem=sd(value)/sqrt(length(value)),
                     lower = activity - sem,
                     upper = activity + sem) %>%
    dplyr::filter(roi_file=="rCA1_body" | roi_file=="rCA2_3_DG_body",
                  condition %in% c("brownRHitsxFHits_Miss", "grayRHitsxFHits_Miss"),
                  hemi == "ashs_left") %>%
    # re-label to pretty up plotting
    dplyr::mutate(roi_lbl = car::recode(roi_file, "'rCA1_body' = 'left CA1'; 'rCA2_3_DG_body' = 'left CA23DG'"),
                  condition_lbl = car::recode(condition, "'brownRHitsxFHits_Miss' = 'Brown House'; 'grayRHitsxFHits_Miss' = 'Gray House'")) %>%
    ggplot2::ggplot(ggplot2::aes(condition_lbl, activity, fill=condition_lbl)) +
    ggplot2::geom_bar(width=0.7,position=ggplot2::position_dodge(0.9),stat="identity") +
    ggplot2::facet_grid(.~roi_lbl) +
    ggplot2::ylab("Mean Univariate Activity") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lower,ymax=upper),position=ggplot2::position_dodge(width=0.9),color="black",width=0.15) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, color = "black"), axis.title.x = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(size = 17),
                   axis.text.y = ggplot2::element_text(size = 10), axis.title.y = ggplot2::element_text(size = 20),
                   strip.text.y = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_blank(),
                   plot.title = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = "none")

  if(SAVE_GRAPHS_FLAG == 1){
    ggplot2::ggsave(file = file.path(graph_fpath_out,
                                     "RHits_FHitsANDMisses_brownVSgray_left-hemi.pdf"),
                    width=8, height=6)
  }

  # --- statiscally compare brown vs gray house activity in the subfields ---
  cons %>%
    dplyr::filter(roi_file == "rCA1_body",
                  contrast %in% c("brownRHitsxFHits_Miss", "grayRHitsxFHits_Miss"),
                  hemi == "ashs_left") %>%
    ez::ezANOVA(.,
                dv=.(activity),
                wid=.(subject),
                within=.(contrast),
                detailed=TRUE)

  cons %>%
    dplyr::filter(roi_file == "rCA2_3_DG_body",
                  contrast %in% c("brownRHitsxFHits_Miss", "grayRHitsxFHits_Miss"),
                  hemi == "ashs_left") %>%
    ez::ezANOVA(.,
                dv=.(activity),
                wid=.(subject),
                within=.(contrast),
                detailed=TRUE)

}
