#' ---
#'  title: Load ABCDCon data
#'  author: Halle R. Dimsdale-Zucker
#'  output:
#'    html_document:
#'      toc: true
#'      number_sections: true
#'      theme: spacelab
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
raw_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$raw_behavioral))
analyzed_behavioral_dir <- paste0(project_dir,halle::ensure_trailing_slash(config$directories$analyzed_behavioral))

#' ## Setup other variables
subjects <- c(1:30)
exclude_subjects <- c(3:5,15:17)
subjects <- subjects[!is.element(subjects, exclude_subjects)]
SAVE_ONSETS_FLAG <- 1

#' # Load data
for(sub in subjects) {

  subNum <- halle::format_three_digit_subject_id(sub)
  sub_objRec <- read.csv(objectRecogFilePath(subNum,raw_behavioral_dir))
  head(sub_objRec)

  # loop across runs and insert column w/ trial numbers
  for(irun in unique(sub_objRec$BlockID)) {

    # if this is the first time through the for loop,
    # initialize sub_objRec_ordered as an empty dataframe
    if(irun==1){
      sub_objRec_ordered <- data.frame()
    }

    sub_objRec_ordered <-
      sub_objRec %>%
      # sort by run and onset w/in that block
      dplyr::arrange(BlockID,Onset) %>%
      # subset current run
      dplyr::filter(BlockID == irun) %>%
      # assign trial numbers for current run
      dplyr::add_rownames() %>% #NB, this is now deprecated in `dplyr`
      # append new rows to sub_objRec_ordered
      dplyr::bind_rows(sub_objRec_ordered,.)
  } #for (irun in

  head(sub_objRec_ordered)

  # overwrite sub_objRec w/ variable that includes trial numbers
  sub_objRec <- sub_objRec_ordered

  # rename "rownames" column
  sub_objRec <-
    sub_objRec %>%
    dplyr::rename(objRec_trialNum = rowname)

  head(sub_objRec)

  # read in location recognition data for current subject
  sub_locRec <- read.csv(locationRecogFilePath(subNum,raw_behavioral_dir))

  # merge objectrec and locrec data for sub (do in multiple steps bc can only merge 2 data frames at a time)
  sub_Rdata <- merge(sub_objRec,sub_locRec, by = c("SubjectID","ObjectIDinMovie","ObjectID","VideoID"),suffixes = c('.objrec','.locrec'),all=TRUE)

  # add to group dataframe
  # if currently working on the first subject, create "data" variable
  # TODO: do not use `data` as a variable name
  if(sub==subjects[1]){
    data <- sub_Rdata
  # if not first subject, bind current subject with "data"
  } else {
    data <- rbind(data,sub_Rdata)
  }
} #for (sub in subjects

#' ## Take a peak at the group data frame
head(data)

#' # Tidy data
# some notes on abbrevs used in variable names:
# HC = house correct
# RC = room correct
# PC = position correct

#' ## fix study locations (CorrectLocationResp) that were incorrect in an earlier version of the StimList_ByObject....txt file
data$CorrectLocationResp[data$ObjectID == "Object205"] <- 15
data$CorrectLocationResp[data$ObjectID == "Object231"] <- 25

#' ## rename, change variables to factors, make new variables, etc.
remove_object_func <- function(s) unlist(strsplit(s, split='t', fixed=TRUE))[2] #eliminate "Object" from "Object###" so just have numeric code

data <-
  data %>%
  # rename unclear variables
  dplyr::rename(studyLocationRCPC = CorrectLocationResp) %>%
  # change variables to be factors
  # based on: http://stackoverflow.com/questions/27668266/dplyr-change-many-data-types
  dplyr::mutate_each(funs(factor),one_of("studyLocationRCPC","CorrectMemResp")) %>%
  # make new variables
  # car::recode tips from: http://rprogramming.net/recode-data-in-r/
  dplyr::mutate(CorrectMemResp_recoded = car::recode(CorrectMemResp,"0='New'; 1='Studied'"),
                LocationResp.byHouse = car::recode(LocationResp, "'1!'='BrownRm1';'2@'='BrownRm2';'3#'='GrayRm1';'4$'='GrayRm2'"),
                LocationResp.byHouse_noRoom = car::recode(LocationResp, "'1!'='Brown';'2@'='Brown';'3#'='Gray';'4$'='Gray'"),
                subj_num = SubjectID,
                # returns 1 when SubjectID odd, 0 when even
                # (just like matlab mod function, which is how counterbalance order was determined)
                CBorder = factor(subj_num %% 2),
                # TODO: seems like there should be a way to do this with `stringr`, but doing this for now
                # figure out which house an item was studied in
                # first do Brown house trials
                studyLocationHC = car::recode(VideoID, "c('Brown1','Brown2','Brown3','Brown4','Brown5','Brown6','Brown7','Brown8','Brown9','Brown10','Brown11','Brown12') = 'Brown'"),
                # now go Gray house trials
                studyLocationHC = car::recode(studyLocationHC, "c('Gray1','Gray2','Gray3','Gray4','Gray5','Gray6','Gray7','Gray8','Gray9','Gray10','Gray11','Gray12') = 'Gray'"),
                # figure out which room an item was studied in
                studyLocationRC = car::recode(studyLocationRCPC, "11:18 = 'Rm1'; 21:28 = 'Rm2'"),
                # figure out which house *and* room an item was studied in
                studyLocationHCRC = paste0(studyLocationHC,studyLocationRC),
                # figure out the house, room, and position an item was studied in
                studyLocationHCRCPC = paste0(studyLocationHC, "Rm", studyLocationRCPC),
                # eliminate "Object" from "Object###" so just have numeric code
                ObjectID_numbers = sapply(as.character(ObjectID),remove_object_func)) %>%
  # and make some more factors
  # (these are things that cannot be set as factors earlier b/c are used to make other variables or hadn't been created by then)
  dplyr::mutate_each(funs(factor),one_of("SubjectID",
                                         "CorrectMemResp",
                                         "studyLocationHC",
                                         "studyLocationRC",
                                         "studyLocationRCPC",
                                         "studyLocationHCRCPC")) %>%
  # set things as factors that are picky about their levels
  dplyr::mutate(studyLocationHCRC = factor(studyLocationHCRC, levels=c("BrownRm1","BrownRm2","GrayRm1","GrayRm2")),
                LocationResp.byHouse = factor(LocationResp.byHouse,levels=c("BrownRm1","BrownRm2","GrayRm1","GrayRm2",exclude=c("5%","xx"))),
                LocationResp.byHouse_noRoom = factor(LocationResp.byHouse_noRoom, levels = c("Brown", "Gray", exclude = c("5%","xx")))) %>%
  # remove duplicate variables (from revalue-ing in mutate)
  dplyr::select(-CorrectMemResp, -ObjectID) %>%
  # rename variables (from revalue-ing)
  dplyr::rename(CorrectMemResp = CorrectMemResp_recoded,
                ObjectID = ObjectID_numbers)

#' ## Recode responses by CB
data1 <- subset(data,CBorder==1)
data2 <- subset(data,CBorder==0)
data1$MemResp <- car::recode(data1$MemResp,c("'1!'='R'; '2@'='F'; '3#'='N'"))
data2$MemResp <- car::recode(data2$MemResp,c("'1!'='N'; '2@'='F';'3#'='R'"))
data <- rbind(data1,data2)

#' ### figure out name labels, by CB order
data$LocationResp.byName <- factor(NA,levels=c("AlexRm1","AlexRm2","JamieRm1","JamieRm2")) #create a new column (as a placeholder) filled w NAs

# for ODD subjects: name1 = Alex, name2= Jamie (from context enc), for EVEN subjects: name1 = Jamie, name2= Alex (from context enc)
for(sub in subjects) { #TODO: Iterate across CB order instead of looping across subjects

    currdata <- subset(data,data$subj_num == sub)
    House1Idx <- which(!is.na(currdata$House1))

    if (currdata$CBorder[1] == 1){
      if (grepl(currdata$House1[House1Idx[1]],"Model_25_Brown_HiddenLayers_BirdsEyeAtEnd")){
        currdata$LocationResp.byName <- car::recode(currdata$LocationResp,c("'1!'='AlexRm1';'2@'='AlexRm2';'3#'='JamieRm1';'4$'='JamieRm2'"))
      } else if (grepl(currdata$House1[House1Idx[1]],"Model_25_Gray_HiddenLayers_BirdsEyeAtEnd")){
        currdata$LocationResp.byName <- car::recode(currdata$LocationResp,c("'1!'='JamieRm1';'2@'='JamieRm2';'3#'='AlexRm1';'4$'='AlexRm2'"))
      }
    } else if (currdata$CBorder[1] == 0){
      if (grepl(currdata$House1[House1Idx[1]],"Model_25_Brown_HiddenLayers_BirdsEyeAtEnd")){
        currdata$LocationResp.byName <- car::recode(currdata$LocationResp,c("'1!'='JamieRm1';'2@'='JamieRm2';'3#'='AlexRm1';'4$'='AlexRm2'"))
      } else if (grepl(currdata$House1[House1Idx[1]],"Model_25_Gray_HiddenLayers_BirdsEyeAtEnd")){
        currdata$LocationResp.byName <- car::recode(currdata$LocationResp,c("'1!'='AlexRm1';'2@'='AlexRm2';'3#'='JamieRm1';'4$'='JamieRm2'"))
      }
    }

    data[which(data$subj_num == sub),] <- currdata
    #NB: This will return a warning about "invalid factor level, NA generated" (and it's okay because not every subject has every combination from if loops above--that's how CB works!)
}

#' ## Make char responses factors w/ levels
# These cannot be done earlier because some of these values get changed in the above for loop based on counterbalance order
data <-
  data %>%
  dplyr::mutate(LocationResp.byName = factor(LocationResp.byName, levels = c("AlexRm1","AlexRm2","JamieRm1","JamieRm2")),
                MemResp = factor(MemResp, levels = c("R","F","N", exclude = c("xx","4$"))))

#' ## Again, take a peak at the data now that it's been tidyed
head(data)

#' # Create columns that score location responses (HC= house correct, HCRC = house correct, room correct)
data$LocationResp.HC[data$studyLocationHC == "Brown" & (data$LocationResp.byHouse=="BrownRm1"|data$LocationResp.byHouse=="BrownRm2")|
                            data$studyLocationHC == "Gray" & (data$LocationResp.byHouse=="GrayRm1"|data$LocationResp.byHouse=="GrayRm2")] <- 1
data$LocationResp.HC[is.na(data$LocationResp.HC)] <- 0

data$LocationResp.HCRC[(data$studyLocationHCRC == "BrownRm1" & data$LocationResp.byHouse=="BrownRm1")|
                            (data$studyLocationHCRC == "BrownRm2" & data$LocationResp.byHouse=="BrownRm2")|
                            (data$studyLocationHCRC == "GrayRm1" & data$LocationResp.byHouse=="GrayRm1")|
                            (data$studyLocationHCRC == "GrayRm2" & data$LocationResp.byHouse=="GrayRm2")] <- 1
data$LocationResp.HCRC[is.na(data$LocationResp.HCRC)] <- 0

#' ## And again, look at the data now that location responses have been scored
head(data)

#' # Create onsets files
#' ## Create regressor numbers
# Based on combination of memory response (data$MemResp), correct memory response (data$CorrectMemResp),...
# location response (data$LocationResp.byHouse), and correct study location (data$studyLocationHCRC)

#' ### memResp hits, locationRespHits (within house, wihtin room (HCRC))
data$trialType[data$LocationResp.byHouse == "BrownRm1" & data$studyLocationHCRC == "BrownRm1" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 1
data$trialType[data$LocationResp.byHouse == "BrownRm2" & data$studyLocationHCRC == "BrownRm2" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 2
data$trialType[data$LocationResp.byHouse == "GrayRm1" & data$studyLocationHCRC == "GrayRm1" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 3
data$trialType[data$LocationResp.byHouse == "GrayRm2" & data$studyLocationHCRC == "GrayRm2" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 4
data$trialType[data$LocationResp.byHouse == "BrownRm1" & data$studyLocationHCRC == "BrownRm1" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 5
data$trialType[data$LocationResp.byHouse == "BrownRm2" & data$studyLocationHCRC == "BrownRm2" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 6
data$trialType[data$LocationResp.byHouse == "GrayRm1" & data$studyLocationHCRC == "GrayRm1" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 7
data$trialType[data$LocationResp.byHouse == "GrayRm2" & data$studyLocationHCRC == "GrayRm2" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 8

#' ### memResp hits, locationResp within house misses (HCRI)
data$trialType[data$LocationResp.byHouse == "BrownRm1" & data$studyLocationHCRC == "BrownRm2" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 9
data$trialType[data$LocationResp.byHouse == "BrownRm2" & data$studyLocationHCRC == "BrownRm1" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 10
data$trialType[data$LocationResp.byHouse == "GrayRm1" & data$studyLocationHCRC == "GrayRm2" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 11
data$trialType[data$LocationResp.byHouse == "GrayRm2" & data$studyLocationHCRC == "GrayRm1" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 12
data$trialType[data$LocationResp.byHouse == "BrownRm1" & data$studyLocationHCRC == "BrownRm2" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 13
data$trialType[data$LocationResp.byHouse == "BrownRm2" & data$studyLocationHCRC == "BrownRm1" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 14
data$trialType[data$LocationResp.byHouse == "GrayRm1" & data$studyLocationHCRC == "GrayRm2" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 15
data$trialType[data$LocationResp.byHouse == "GrayRm2" & data$studyLocationHCRC == "GrayRm1" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 16

#' ### memResp hits, locationResp (HI)
data$trialType[data$LocationResp.byHouse == "BrownRm1" & data$studyLocationHC == "Gray" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 17
data$trialType[data$LocationResp.byHouse == "BrownRm2" & data$studyLocationHC == "Gray" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 17
data$trialType[data$LocationResp.byHouse == "GrayRm1" & data$studyLocationHC == "Brown" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 18
data$trialType[data$LocationResp.byHouse == "GrayRm2" & data$studyLocationHC == "Brown" & data$MemResp == "R" & data$CorrectMemResp == "Studied"] <- 18
data$trialType[data$LocationResp.byHouse == "BrownRm1" & data$studyLocationHC == "Gray" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 19
data$trialType[data$LocationResp.byHouse == "BrownRm2" & data$studyLocationHC == "Gray" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 19
data$trialType[data$LocationResp.byHouse == "GrayRm1" & data$studyLocationHC == "Brown" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 20
data$trialType[data$LocationResp.byHouse == "GrayRm2" & data$studyLocationHC == "Brown" & data$MemResp == "F" & data$CorrectMemResp == "Studied"] <- 20

#' ### memResp misses
data$trialType[data$MemResp == "N" & data$CorrectMemResp == "Studied"] <- 21

#' ### memResp CR
data$trialType[data$MemResp == "N" & data$CorrectMemResp == "New"] <- 22

#' ### memResp FA
data$trialType[data$MemResp == "R" & data$CorrectMemResp == "New"] <- 23
data$trialType[data$MemResp == "F" & data$CorrectMemResp == "New"] <- 24

#' ### trials of no interest (ie, not included)
data$trialType[data$LocationResp == "xx"] <- 99
data$trialType[data$LocationResp == "5%"] <- 99

#' ### no memResp
data$trialType[data$MemResp == "xx"] <- 99
data$trialType[data$MemResp == "4$"] <- 99

#' ### key
data$trialTypeKey <- car::recode(factor(data$trialType),c("1='BrownRm1_RHit_HCRC';2='BrownRm2_RHit_HCRC';3='GrayRm1_RHit_HCRC';4='GrayRm2_RHit_HCRC';
                                                            5='BrownRm1_FHit_HCRC';6='BrownRm2_FHit_HCRC';7='GrayRm1_FHit_HCRC';8='GrayRm2_FHit_HCRC';
                                                            9='BrownRm2_RHit_HCRI';10='BrownRm1_RHit_HCRI';11='GrayRm2_RHit_HCRI';12='GrayRm1_RHit_HCRI';
                                                            13='BrownRm2_FHit_HCRI';14='BrownRm1_FHit_HCRI';15='GrayRm2_FHit_HCRI';16='GrayRm1_FHit_HCRI';
                                                            17='Gray_RHit_HI';18='Brown_RHit_HI';
                                                            19='Gray_FHit_HI';20='Brown_FHit_HI';
                                                            21='Miss';22='CR';23='FA_R';24='FA_F';
                                                            99='ExcludeTrial'"))

levels(data$trialTypeKey)

#' ## Switch to "o" dataframe
# This way, onsets information is separate from what gets saved out in group file
o <- NULL
o$subj_num <- data$subj_num
o$trialType <- data$trialType
o$trialTypeKey <- data$trialTypeKey
o$Onset <- data$Onset
o$BlockID <- data$BlockID
o$ObjectID <- data$ObjectID
o$VideoID <- data$VideoID
o$duration <- 0
o$JitterLength <- data$JitterLength
o$objRec_trialNum <- data$objRec_trialNum
o <- data.frame(o)

#' ## Sanity check model
#' ### Revalue regressors for sanity check model
# This should include house1_Rhits, house2_Rhits, house1_Fhits&misses, house2_Fhits&misses, CR, FA, NoResp
# Regressor numbers from above will be collapsed together to fit into these current regressors of interest
o$sanityCheck_model <- car::recode(factor(data$trialTypeKey),c("'BrownRm1_RHit_HCRC'=1;'BrownRm2_RHit_HCRC'=1;'BrownRm2_RHit_HCRI'=1;'BrownRm1_RHit_HCRI'=1;'Brown_RHit_HI'=1;
                                                               'GrayRm1_RHit_HCRC'=2;'GrayRm2_RHit_HCRC'=2;'GrayRm2_RHit_HCRI'=2;'GrayRm1_RHit_HCRI'=2;'Gray_RHit_HI'=2;
                                                               'BrownRm1_FHit_HCRC'=3;'BrownRm2_FHit_HCRC'=3;'BrownRm1_FHit_HCRI'=3;'BrownRm2_FHit_HCRI'=3;'Brown_FHit_HI'=3;
                                                               'GrayRm1_FHit_HCRC'=4;'GrayRm2_FHit_HCRC'=4;'GrayRm1_FHit_HCRI'=4;'GrayRm2_FHit_HCRI'=4;'Gray_FHit_HI' =4;
                                                               'CR'=5;
                                                               'FA_R'=6;'FA_F'=6;
                                                               'ExcludeTrial'=99"))
o$sanityCheck_model[data$trialTypeKey=="Miss" & data$studyLocationHC=="Brown"] <- 3
o$sanityCheck_model[data$trialTypeKey=="Miss" & data$studyLocationHC=="Gray"] <- 4
# eliminate "Miss" from being counted as a level
o$sanityCheck_model <- droplevels(o$sanityCheck_model)

#' ### Revalue "key" for regressors in sanity check model
o$sanityCheck_key <- car::recode(factor(data$trialTypeKey),c("'BrownRm1_RHit_HCRC'='Brown_RHit';'BrownRm2_RHit_HCRC'='Brown_RHit';'BrownRm2_RHit_HCRI'='Brown_RHit';'BrownRm1_RHit_HCRI'='Brown_RHit';'Brown_RHit_HI'='Brown_RHit';
                                                             'GrayRm1_RHit_HCRC'='Gray_RHit';'GrayRm2_RHit_HCRC'='Gray_RHit';'GrayRm2_RHit_HCRI'='Gray_RHit';'GrayRm1_RHit_HCRI'='Gray_RHit';'Gray_RHit_HI'='Gray_RHit';
                                                             'BrownRm1_FHit_HCRC'='Brown_FHitsANDMisses';'BrownRm2_FHit_HCRC'='Brown_FHitsANDMisses';'BrownRm1_FHit_HCRI'='Brown_FHitsANDMisses';'BrownRm2_FHit_HCRI'='Brown_FHitsANDMisses';'Brown_FHit_HI'='Brown_FHitsANDMisses';
                                                             'GrayRm1_FHit_HCRC'='Gray_FHitsANDMisses';'GrayRm2_FHit_HCRC'='Gray_FHitsANDMisses';'GrayRm1_FHit_HCRI'='Gray_FHitsANDMisses';'GrayRm2_FHit_HCRI'='Gray_FHitsANDMisses';'Gray_FHit_HI' ='Gray_FHitsANDMisses';
                                                             'FA_R'='FA';'FA_F'='FA'"))
o$sanityCheck_key[data$trialTypeKey=="Miss" & data$studyLocationHC=="Brown"] <- "Brown_FHitsANDMisses"
o$sanityCheck_key[data$trialTypeKey=="Miss" & data$studyLocationHC=="Gray"] <- "Gray_FHitsANDMisses"
# eliminate "Miss" from being counted as a level
o$sanityCheck_key <- droplevels(o$sanityCheck_key)

#' ## Model just by memory, ignoring all contexts
# This should include RHits, Fhits, Misses, CR, FA, and ExcludeTrials
o$byMemory_model <- car::recode(factor(data$trialTypeKey),c("'BrownRm1_RHit_HCRC'=1;'BrownRm2_RHit_HCRC'=1;'BrownRm2_RHit_HCRI'=1;'BrownRm1_RHit_HCRI'=1;'Brown_RHit_HI'=1;
                                                               'GrayRm1_RHit_HCRC'=1;'GrayRm2_RHit_HCRC'=1;'GrayRm2_RHit_HCRI'=1;'GrayRm1_RHit_HCRI'=1;'Gray_RHit_HI'=1;
                                                               'BrownRm1_FHit_HCRC'=2;'BrownRm2_FHit_HCRC'=2;'BrownRm1_FHit_HCRI'=2;'BrownRm2_FHit_HCRI'=2;'Brown_FHit_HI'=2;
                                                               'GrayRm1_FHit_HCRC'=2;'GrayRm2_FHit_HCRC'=2;'GrayRm1_FHit_HCRI'=2;'GrayRm2_FHit_HCRI'=2;'Gray_FHit_HI' =2;
                                                               'Miss'=3;
                                                               'CR'=4;
                                                               'FA_R'=5;'FA_F'=5;
                                                               'ExcludeTrial'=99"))
#' ### Revalue "key" for regressors in memory model
# no need to recode misses or excluded trials b/c they already have the correct labels
o$byMemory_key <- car::recode(factor(data$trialTypeKey),c("'BrownRm1_RHit_HCRC'='RHit';'BrownRm2_RHit_HCRC'='RHit';'BrownRm2_RHit_HCRI'='RHit';'BrownRm1_RHit_HCRI'='RHit';'RHit_HI'='RHit';'Brown_RHit_HI'='RHit';
                                                             'GrayRm1_RHit_HCRC'='RHit';'GrayRm2_RHit_HCRC'='RHit';'GrayRm2_RHit_HCRI'='RHit';'GrayRm1_RHit_HCRI'='RHit';'RHit_HI'='RHit';'Gray_RHit_HI'='RHit';
                                                             'BrownRm1_FHit_HCRC'='FHit';'BrownRm2_FHit_HCRC'='FHit';'BrownRm1_FHit_HCRI'='FHit';'BrownRm2_FHit_HCRI'='FHit';'Brown_FHit_HI'='FHit';
                                                             'GrayRm1_FHit_HCRC'='FHit';'GrayRm2_FHit_HCRC'='FHit';'GrayRm1_FHit_HCRI'='FHit';'GrayRm2_FHit_HCRI'='FHit';'Gray_FHit_HI' ='FHit';
                                                             'FA_R'='FA';'FA_F'='FA'"))
levels(o$byMemory_key)

#' ## Model just by memory, combining Fhits & Misses, ignoring all contexts
# This should include RHits, FhitsANDMisses, CR, FA, and ExcludeTrials
o$byMemoryFHitMiss_model <- car::recode(factor(data$trialTypeKey),c("'BrownRm1_RHit_HCRC'=1;'BrownRm2_RHit_HCRC'=1;'BrownRm2_RHit_HCRI'=1;'BrownRm1_RHit_HCRI'=1;'Brown_RHit_HI'=1;
                                                                    'GrayRm1_RHit_HCRC'=1;'GrayRm2_RHit_HCRC'=1;'GrayRm2_RHit_HCRI'=1;'GrayRm1_RHit_HCRI'=1;'Gray_RHit_HI'=1;
                                                                    'BrownRm1_FHit_HCRC'=2;'BrownRm2_FHit_HCRC'=2;'BrownRm1_FHit_HCRI'=2;'BrownRm2_FHit_HCRI'=2;'Brown_FHit_HI'=2;
                                                                    'GrayRm1_FHit_HCRC'=2;'GrayRm2_FHit_HCRC'=2;'GrayRm1_FHit_HCRI'=2;'GrayRm2_FHit_HCRI'=2;'Gray_FHit_HI' =2;
                                                                    'Miss'=2;
                                                                    'CR'=3;
                                                                    'FA_R'=4;'FA_F'=4;
                                                                    'ExcludeTrial'=99"))
#' ### Revalue "key" for regressors in memory model
# no need to recode misses or excluded trials b/c they already have the correct labels
o$byMemoryFHitMiss_key <- car::recode(factor(data$trialTypeKey),c("'BrownRm1_RHit_HCRC'='RHit';'BrownRm2_RHit_HCRC'='RHit';'BrownRm2_RHit_HCRI'='RHit';'BrownRm1_RHit_HCRI'='RHit';'RHit_HI'='RHit';'Brown_RHit_HI'='RHit';
                                                                  'GrayRm1_RHit_HCRC'='RHit';'GrayRm2_RHit_HCRC'='RHit';'GrayRm2_RHit_HCRI'='RHit';'GrayRm1_RHit_HCRI'='RHit';'RHit_HI'='RHit';'Gray_RHit_HI'='RHit';
                                                                  'BrownRm1_FHit_HCRC'='FHitsANDMisses';'BrownRm2_FHit_HCRC'='FHitsANDMisses';'BrownRm1_FHit_HCRI'='FHitsANDMisses';'BrownRm2_FHit_HCRI'='FHitsANDMisses';'Brown_FHit_HI'='FHitsANDMisses';
                                                                  'GrayRm1_FHit_HCRC'='FHitsANDMisses';'GrayRm2_FHit_HCRC'='FHitsANDMisses';'GrayRm1_FHit_HCRI'='FHitsANDMisses';'GrayRm2_FHit_HCRI'='FHitsANDMisses';'Gray_FHit_HI' ='FHitsANDMisses';
                                                                  'FA_R'='FA';'FA_F'='FA'"))
o$byMemoryFHitMiss_key[data$trialTypeKey=="Miss"] <- "FHitsANDMisses"
# eliminate "Miss" from being counted as a level
o$byMemoryFHitMiss_key <- droplevels(o$byMemoryFHitMiss_key)
levels(o$byMemoryFHitMiss_key)

#' # Save out
#' ## Onsets
# fix formatting so plays nicely when save out
o <- as.data.frame(o)

if (SAVE_ONSETS_FLAG==1){
  for(sub in subjects){
    if(is.element(sub, c(3,15))){
      currdata <- subset(o,o$subj_num == sub)
      # eliminate 'NA' rows
      currdata <- currdata[complete.cases(currdata),]
    } else {
      currdata <- subset(o,o$subj_num == sub)
    }

    subNum <- halle::format_three_digit_subject_id(sub)
    raw_behavioral_dir <- halle::ensure_trailing_slash(raw_behavioral_dir)
    subjDir <- paste0(raw_behavioral_dir, subNum)
    subjDir <- halle::ensure_trailing_slash(subjDir)
    dir.create(subjDir, recursive = TRUE)
    write.csv(currdata,file=paste0(subjDir,'onsets_file.csv'),row.names=FALSE,quote=FALSE)
  }
}

#'## Group files
# save out data to .RData and .csv files
save(data,file=paste0(raw_behavioral_dir,'group_data.RData'))
write.csv(data,file=paste0(raw_behavioral_dir,'group_data.csv'))

paste(cat("Current number of subjects included: "),length(unique(data$subj_num)),sep="")

