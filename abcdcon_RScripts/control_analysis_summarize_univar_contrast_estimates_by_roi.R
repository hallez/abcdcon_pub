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

#' ## Set paths as variables
analyzed_mri_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_mri))
raw_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$raw_behavioral))
analyzed_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_behavioral))
contrast_est_path <- paste0(analyzed_mri_dir,halle::ensure_trailing_slash('univariate_sanityCheck'))
contrast_est_fname <- 'cons_body_ROIs.csv' # this is created by `control_analysis_univariate_contrast_estimates_by_roi.m`
contrast_est_fpath <- paste0(contrast_est_path,contrast_est_fname)

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
    dplyr::filter(roi_file=="rCA1_body" | roi_file=="rCA2_3_DG_body" | roi_file=="rwhole_hippo",
                  condition %in% c("brownRHitsxFHits_Miss", "grayRHitsxFHits_Miss"),
                  hemi == "ashs_left") %>%
    ggplot2::ggplot(ggplot2::aes(roi_file, activity, fill=roi_file)) +
    ggplot2::geom_bar(width=0.7,position=ggplot2::position_dodge(0.9),stat="identity") +
    ggplot2::facet_grid(.~condition) +
    ggplot2::ylab("mean activity") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), strip.text.x = ggplot2::element_text(size = 15),
                   axis.text.y = ggplot2::element_text(size = 20), axis.title.y=ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size=20)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lower,ymax=upper),position=ggplot2::position_dodge(width=0.9),color="black",width=0.25)
}
