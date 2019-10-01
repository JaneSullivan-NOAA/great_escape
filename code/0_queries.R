# Query and clean data - processing script for all incoming data
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-10-01

# All final columns formatted as follows unless otherwise specified:
# length: fork length, cm
# weight: kg
# depth: average depth, meters
# date: ISO 8601 YYYY/MM/DD
# characters & factors - first letter capitilized (e.g. 'Sex'), otherwise lowercase
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

# Queries -----

# Experimental code for escape ring study = 66 
query <- 
  paste0(
" select  year, effort_no, 
          length_millimeters, weight_kilograms, girth_millimeters,
          pot_treatment_code, pot_treatment
                   
  from    out_g_bio_effort_age_sex_size

  where   species_code = '710' and
          project_code in ('66') and
          year = ", YEAR)

dbGetQuery(zprod_channel, query) -> pot_bio

write_csv(pot_bio, paste0("data/raw/pot_bio_", YEAR, ".csv"))

read_csv(paste0("data/raw/pot_bio_", YEAR, ".csv"), guess_max = 50000) %>% 
  mutate(Treatment = derivedFactor("Control" = POT_TREATMENT_CODE == "99",
                                   "4.00 in" = POT_TREATMENT_CODE == "00",
                                   "3.75 in" = POT_TREATMENT_CODE == "01",
                                   "3.50 in" = POT_TREATMENT_CODE == "02",
                                   .default = NA),
         length_cm = LENGTH_MILLIMETERS / 10) %>% 
  select(year = YEAR, effort_no = EFFORT_NO, Treatment, length_cm, weight_kg = WEIGHT_KILOGRAMS,
         girth_mm = GIRTH_MILLIMETERS) %>% 
  write_csv(paste0("data/pot_bio_", YEAR, ".csv"))
