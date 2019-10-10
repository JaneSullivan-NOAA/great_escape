# Exploratory data analysis
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-10-01

# Set up ----

YEAR <- 2019 # study year(s)
source("code/helper.r")

# Bio data
bio <- read_csv(paste0("data/pot_bio_", YEAR, ".csv")) %>% 
  filter(!is.na(length)) %>% 
  mutate(Treatment = derivedFactor("Control" = Treatment == "99",
                                   "3.50 in" = Treatment == "02",
                                   "3.75 in" = Treatment == "01",
                                   "4.00 in" = Treatment == "00",
                                   .default = NA,
                                   .ordered = TRUE))
# Email sent to A. Baldwin 2019-10-01 about the NAs. See GitHub issue #1
bio <- bio %>% filter(!is.na(Treatment))

# Effort data 
effort <- read_csv(paste0("data/pot_effort_", YEAR, ".csv"))

install.packages("glmmTMB")
library(glmmTMB)

mod <- 

# EDA ----



# The control pot captured 

