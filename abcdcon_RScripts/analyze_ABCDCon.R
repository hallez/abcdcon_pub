#' ---
#'  title: Analyze ABCDCon
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
dropbox_dir <- paste0(halle::ensure_trailing_slash(config$directories$dropbox_abcdcon))

#' ## Setup other variables
SAVE_GRAPHS_FLAG <- 1
TESTING_FLAG <- 0

#' # Load group data
load(paste0(raw_behavioral_dir,'group_data.RData'))
str(data)

#' # Only include subjects who are in final MRI analyses
# ideally, wouldn't overwrite variable w/ the same name
# (also, ideally `data` wouldn't be the variable name)
data <- data %>%
  dplyr::filter(subj_num %in% c(1:2, 6:14, 18:21, 23:30))

#' # Object recognition
#' ## Rhits, by subject
#' TODO: figure out how to show summary output for all subjects
data_objrec_summarized <-
  data %>%
  # pull out only relevant columns of data
  dplyr::select(subj_num,CorrectMemResp,MemResp) %>%
  # filter so just get R hits
  dplyr::filter(CorrectMemResp=="Studied", MemResp=="R") %>%
  # group by relevant factors
  dplyr::group_by(subj_num,CorrectMemResp,MemResp) %>%
  dplyr::summarise(count = n())

data_objrec_summarized

#' ## Response rates
options(dplyr.width = Inf)
data_objrec_rates <-
  data %>%
  # select only necessary columns
  dplyr::select(subj_num, CorrectMemResp, MemResp) %>%
  # drop weird responses ($, xx, etc.)
  dplyr::filter(MemResp=="R" | MemResp =="F" | MemResp == "N") %>%
  # compute count per response and correct bin
  dplyr::group_by(subj_num, MemResp, CorrectMemResp) %>%
  dplyr::summarize(count = n()) %>%
  # add in zeros for those who have no responses of that type
  # based on http://stackoverflow.com/questions/25811756/summarizing-counts-of-a-factor-with-dplyr
  dplyr::ungroup() %>%
  tidyr::spread(MemResp, count, fill = 0) %>%
  tidyr::gather(MemResp, count, R:N) %>%
  # transform so have one row per subject
  tidyr::unite(memresp_correctresp, MemResp, CorrectMemResp) %>%
  tidyr::spread(memresp_correctresp, count) %>%
  # replace NA values with 0s - this is necessary to compute recollection and familiarity
  dplyr::ungroup() %>%
  dplyr::mutate_each(funs(replace_NA_with_zeros)) %>%
  # compute total number of trials per subject
  dplyr::group_by(subj_num) %>%
  # compute rates
  dplyr::mutate(total_trials = sum(F_New, F_Studied,
                                   N_New, N_Studied,
                                   R_New, R_Studied,
                                   na.rm = TRUE),
                total_old_trials = sum(R_Studied,
                                           F_Studied,
                                           N_Studied,
                                           na.rm = TRUE),
                total_new_trials = sum(R_New,
                                       F_New,
                                       N_New,
                                       na.rm = TRUE),
                r_responses = sum(R_New, R_Studied, na.rm = TRUE),
                f_responses = sum(F_New, F_Studied, na.rm = TRUE),
                n_responses = sum(N_New, N_Studied, na.rm = TRUE),
                r_hit_rate = R_Studied / total_old_trials,
                f_hit_rate = F_Studied / total_old_trials,
                n_miss_rate = N_Studied / total_old_trials,
                r_FA_rate = R_New / total_new_trials,
                f_FA_rate = F_New / total_new_trials,
                n_CR_rate = N_New / total_new_trials)

head(data_objrec_rates)

#' ### Summarize these rates (Supplemental Table 1: Item Recognition)
data_objrec_rates %>%
  ungroup() %>%
  dplyr::select(r_hit_rate, f_hit_rate, n_miss_rate,
                r_FA_rate, f_FA_rate, n_CR_rate) %>%
  dplyr::summarise_each(funs(mean(., na.rm = TRUE),
                      sd(., na.rm = TRUE)))

#' ## Rates by encoding house (Supplemental: Behavioral results)
itemXhouse_respcount <- data %>%
  group_by(subj_num, MemResp, CorrectMemResp, studyLocationHC) %>%
  dplyr::summarise(respcount = length(MemResp))

itemXhouse_rates <- data %>%
  dplyr::group_by(subj_num, CorrectMemResp, studyLocationHC) %>%
  dplyr::summarise(allcount = length(CorrectMemResp)) %>%
  dplyr::left_join(itemXhouse_respcount) %>%
  dplyr::mutate(resprate = respcount / allcount) %>%
  dplyr::group_by(MemResp, CorrectMemResp, studyLocationHC)

itemXhouse_rates %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate),
                   num_subj = length(resprate),
                   sem_resprate = sd_resprate/sqrt(num_subj))

itemXhouse_rates %>%
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::filter(MemResp == "R") %>%
  t.test(resprate ~ studyLocationHC, paired = TRUE, data = .)

# compute Cohen's d
itemXhouse_rates %>%
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::filter(MemResp == "R") %>%
  # remove this so that spread works w/o NA values
  dplyr::select(-respcount) %>%
  tidyr::spread(studyLocationHC, resprate) %>%
  dplyr::mutate(var1 = Brown,
                var2 = Gray) %>%
  dplyr::ungroup() %>%
  dplyr::select(var1, var2) %>%
  halle::compute_cohens_d()

# plot to understand why have significant diff btwn the houses since group rates are similar
itemXhouse_rates %>%
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::filter(MemResp == "R") %>%
  ggplot2::ggplot(., ggplot2::aes(x = studyLocationHC, y = resprate)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(ggplot2::aes(color = as.factor(subj_num))) +
  ggplot2::geom_line(ggplot2::aes(x = studyLocationHC, y = resprate, group = subj_num))

#' # Location recognition
#' ## Test whether response rates varied by house (Supplemental: Behavioral results)
loc_respcount <- data %>%
  # only old items were included in the location source memory test
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::group_by(subj_num, LocationResp.byHouse_noRoom, LocationResp.HC, studyLocationHC) %>%
  dplyr::summarise(respcount = length(LocationResp.byHouse_noRoom))

loc_rates <- data %>%
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::group_by(subj_num, studyLocationHC) %>%
  dplyr::summarise(allcount = length(studyLocationHC)) %>%
  dplyr::left_join(loc_respcount) %>%
  dplyr::mutate(resprate = respcount / allcount) %>%
  dplyr::group_by(LocationResp.byHouse_noRoom, LocationResp.HC)

loc_rates %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate),
                   num_subj = length(resprate),
                   sem_resprate = sd_resprate/sqrt(num_subj))

loc_rates %>%
  dplyr::filter(LocationResp.HC == 1) %>%
  t.test(resprate ~ studyLocationHC, paired = TRUE, data = .)

# compute cohen's d
loc_rates %>%
  dplyr::ungroup() %>%
  dplyr::filter(LocationResp.HC == 1) %>%
  # remove so that spread works w/o NA values
  dplyr::select(-respcount, -LocationResp.byHouse_noRoom) %>%
  tidyr::spread(studyLocationHC, resprate) %>%
  dplyr::mutate(var1 = Brown,
                var2 = Gray) %>%
  dplyr::select(var1, var2) %>%
  halle::compute_cohens_d()

#' ## Compute location rates split by item performance (Supplemental: Behavioral results)
locXitem_resprate <- data %>%
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::group_by(subj_num, LocationResp.byHouse_noRoom, LocationResp.HC, studyLocationHC, MemResp) %>%
  dplyr::summarise(respcount = length(LocationResp.byHouse_noRoom))

locXitem_rates <- data %>%
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::group_by(subj_num, studyLocationHC, MemResp) %>%
  dplyr::summarise(allcount = length(studyLocationHC)) %>%
  dplyr::left_join(locXitem_resprate) %>%
  dplyr::mutate(resprate = respcount / allcount) %>%
  dplyr::group_by(LocationResp.byHouse_noRoom, MemResp)

locXitem_rates %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate),
                   num_subj = length(resprate),
                   sem_resprate = sd_resprate/sqrt(num_subj))

locXitem_rates %>%
  dplyr::filter(LocationResp.byHouse_noRoom %in% c("Brown", "Gray")) %>%
  dplyr::filter(MemResp %in% c("R", "F")) %>%
  dplyr::filter(LocationResp.HC == 1) %>%
  ez::ezANOVA(.,
              dv=.(resprate),
              wid=.(subj_num),
              within=.(MemResp),
              detailed = TRUE)

#' ## Compute source memory w/o breaking down by location (Supplemental Table 1: Spatial Context Source Memory)
loc_collapseHouse_respcount <- data %>%
  # only include trials where subjects made a (sensible) response
  dplyr::filter(LocationResp %in% c("1!", "2@", "3#", "4$")) %>%
  dplyr::group_by(subj_num, LocationResp.HC) %>%
  dplyr::summarise(respcount = length(LocationResp.HC))

loc_collapseHouse_rates <- data %>%
  # only include studied items (since these were the only ones on the location source test)
  dplyr::filter(CorrectMemResp == "Studied") %>%
  dplyr::group_by(subj_num) %>%
  dplyr::summarise(allcount = length(LocationResp.HC)) %>%
  dplyr::left_join(loc_collapseHouse_respcount) %>%
  dplyr::mutate(resprate = respcount / allcount) %>%
  dplyr::group_by(LocationResp.HC)

loc_collapseHouse_rates %>%
  dplyr::summarise(mean_resprate = mean(resprate),
                   sd_resprate = sd(resprate),
                   num_subj = length(resprate),
                   sem_resprate = sd_resprate/sqrt(num_subj))

#' # Compute RTs for object recog (by condition)
recdata <- data %>%
  # drop weird responses ($, xx, etc.)
  dplyr::filter(MemResp=="R" | MemResp =="F" | MemResp == "N")

#' ## Score memory data
oldmask <- which(recdata$CorrectMemResp== "Studied")
newmask <- which(recdata$CorrectMemResp== "New")

Rresp <- which(recdata$MemResp == "R")
Fresp <- which(recdata$MemResp == "F")
Nresp <- which(recdata$MemResp == "N")

recdata$itemScore <- factor(NA,levels=c("Rec","Fam","Miss","R-FA","F-FA","CR"))
recdata$itemScore[intersect(oldmask,Rresp)] <- "Rec"
recdata$itemScore[intersect(oldmask,Fresp)] <- "Fam"
recdata$itemScore[intersect(oldmask,Nresp)] <- "Miss"
recdata$itemScore[intersect(newmask,Rresp)] <- "R-FA"
recdata$itemScore[intersect(newmask,Fresp)] <- "F-FA"
recdata$itemScore[intersect(newmask,Nresp)] <- "CR"

table(recdata$subj_num,recdata$itemScore)

#' ## RTs by item memory
recdata %>%
  dplyr::group_by(itemScore) %>%
  dplyr::summarise(mean_RT = mean(MemRT, na.rm = TRUE))

#' ## RTs by item and source memory
recdata %>%
  dplyr::group_by(itemScore, LocationResp.HC) %>%
  dplyr::summarise(mean_RT = mean(MemRT, na.rm = TRUE))
