#' Create filepath for object recog file
#'
#' @param subNum Subject number. Can use getSubNum.R to format
#' @param rawDataDir Parent directory for raw data

objectRecogFilePath <- function(subNum, rawDataDir) {
  rawDataDir <- halle::ensure_trailing_slash(rawDataDir)
  filePath <- paste0(rawDataDir, subNum,'/ConABCD_objectRecog_',subNum, '.dat')
  return(filePath)
}
