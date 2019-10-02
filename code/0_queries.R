# Query and clean data - processing script for all incoming data
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-10-01

# All final columns formatted as follows unless otherwise specified:
# length: fork length, cm
# girth: mm
# weight: kg
# depth: start, end, average depth - meters
# date: ISO 8601 YYYY/MM/DD
# lat and lon are in decimal degrees

YEAR <- 2019 # study year(s)

# Load ----
source("code/helper.r")

# Oracle connections ----

# Database usernames and passwords (user-specific, ignored)
ora <- read_csv("data/database.csv") 

# Connection strings in T:/Toolbox/TNS/tnsnames.ora

# Zander. Region I database
zprod <- "(DESCRIPTION =
(ADDRESS = (PROTOCOL = TCP)(HOST = 10.209.0.83)(PORT = 1521))
(CONNECT_DATA = (SERVER = DEDICATED)
(SERVICE_NAME = DFGCFR1P.500040564.us1.internal)))"

zprod_channel <- dbConnect(drv = dbDriver('Oracle'), 
                           username = ora$zprod_user, 
                           password = ora$zprod_pw, 
                           dbname = zprod)

# Pot bio data -----

# Experimental project code for escape ring study = 66 
query <- 
  paste0(
" select  year, effort_no, specimen_no,
          length_millimeters, weight_kilograms, girth_millimeters,
          pot_treatment_code, pot_treatment
                   
  from    out_g_bio_effort_age_sex_size

  where   species_code = '710' and
          project_code in ('66') and
          year = ", YEAR)

dbGetQuery(zprod_channel, query) -> pot_bio

write_csv(pot_bio, paste0("data/raw/pot_bio_", YEAR, ".csv"))

read_csv(paste0("data/raw/pot_bio_", YEAR, ".csv"), guess_max = 50000) %>% 
  mutate(treatment = derivedFactor("Control" = POT_TREATMENT_CODE == "99",
                                   "4.00 in" = POT_TREATMENT_CODE == "00",
                                   "3.75 in" = POT_TREATMENT_CODE == "01",
                                   "3.50 in" = POT_TREATMENT_CODE == "02",
                                   .default = NA),
         length = LENGTH_MILLIMETERS / 10) %>% 
  select(year = YEAR, effort_no = EFFORT_NO, specimen_no = SPECIMEN_NO, 
         treatment, length, weight = WEIGHT_KILOGRAMS,
         girth = GIRTH_MILLIMETERS) %>% 
  write_csv(paste0("data/pot_bio_", YEAR, ".csv"))

# Pot effort data ----

# Note this uses project code = 11 (the normal code for the Chatham pot survey)
query <-
  paste0(
" select *

  from out_g_sur_pot

  where species_code = '710' and
        year = ", YEAR)

dbGetQuery(zprod_channel, query) -> pot_effort

write_csv(pot_effort, paste0("data/raw/pot_effort_", YEAR, ".csv"))

read_csv(paste0("data/raw/pot_effort_", YEAR, ".csv"), guess_max = 50000) %>%
  mutate(
    # Define soak time NEEDS_REVIEW
    soak_time = difftime(TIME_FIRST_ANCHOR_ONBOARD, TIME_SECOND_ANCHOR_OVERBOARD, units = "hours"),
    start_depth = START_DEPTH_FATHOMS * 1.8288, # fathoms to meters
    end_depth = END_DEPTH_FATHOMS * 1.8288,
    mean_depth = AVG_DEPTH_FATHOMS * 1.8288) %>% 
  select(year = YEAR, effort_no = EFFORT_NO, stat = G_STAT_AREA,
         start_lat = START_LATITUDE_DECIMAL_DEGREES, start_lon = START_LONGITUDE_DECIMAL_DEGREE,
         end_lat = END_LATITUDE_DECIMAL_DEGREES, end_lon = END_LONGITUDE_DECIMAL_DEGREES,
         soak_time, start_depth, end_depth, mean_depth, n_pots = NUMBER_OF_POTS_RETRIEVED, 
         discard_status_code = DISCARD_STATUS_CODE, discard_status = DISCARD_STATUS,
         n = NUMBERS, comments = POT_COMMENTS) %>% 
  write_csv(paste0("data/pot_effort_", YEAR, ".csv"))