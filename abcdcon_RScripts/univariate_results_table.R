#' ---
#' title: ABCDCon Univariate Results Table
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
project_dir <- "../"
config <- yaml::yaml.load_file(paste0(project_dir,"config.yml"))

#' ## Set paths as variables
analyzed_mri_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_mri))
vendor_dir <- file.path(project_dir,config$directories$vendor)
analysis_name <- "univariate_sanityCheck"
contrast_lbl <- "allRHitsxFHits_Miss_scon0005_N23"

#' # Load in the MNI to AAL function
# remember, this came from: https://github.com/yunshiuan/label4MRI
load(file.path(vendor_dir, "mni2aal.RData"))

#' # Read in ccoordinates from SPM
m <- read.csv(file.path(analyzed_mri_dir, analysis_name, contrast_lbl, "coordinates_p001_k88.csv"), header = TRUE)

#' ## Peek at SPM coordinates
head(m)

#' # Figure out the AAL atlas labels
Result=t(mapply(FUN=mni_to_region_name,x=m$x,y=m$y,z=m$z))

#' ## Split hemi into a separate column
results_fmt <- Result %>%
  as.data.frame() %>%
  tidyr::separate(region, into = c("region_name", "hemi"), sep = "_")

#' ## Print out labels
knitr::kable(results_fmt)

#' # Save out
write.csv(results_fmt,
          file.path(config$directories$dropbox_abcdcon, "writeups", "figures", "univariate_allRHitsxFHits_Misses_p001_k88_aal_lbls.csv"))
