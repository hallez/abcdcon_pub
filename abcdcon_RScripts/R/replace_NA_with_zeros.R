#' Function to replace NAs with 0s
#'
#' based on: http://stackoverflow.com/questions/26081672/change-variable-values-by-groups-using-dplyr
#' @param x Input data that may have NAs

replace_NA_with_zeros <- function(x) {
  ifelse(is.na(x), 0, x)
}
