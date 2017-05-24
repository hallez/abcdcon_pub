#' Create filepath for location recog file
#'
#' @param subNum Subject number. Can use getSubNum.R to format
#' @param rawDataDir Parent directory for raw data

locationRecogFilePath <- function(subNum, rawDataDir) {
  rawDataDir <- halle::ensure_trailing_slash(rawDataDir)
  filePath <- paste0(rawDataDir, subNum,'/ConABCD_locationRecog_',subNum, '.dat')
  return(filePath)
}
