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


# EDA ----

# Length frequencies ----
bio %>% 
  filter(!is.na(Treatment) & between(length, 40, 90)) %>%
  droplevels() %>%
  ggplot(aes(x = length, colour = Treatment, size = Treatment, linetype = Treatment)) + #, fill = treatment)) +
  geom_freqpoly() +
  scale_colour_manual(values = c("grey90", "grey70", "grey40", "black")) +
  scale_size_manual(values = c(1.3, 0.7, 0.7, 0.7)) +
  scale_linetype_manual(values = c(1, 1, 2, 3)) +
  # stat_density_ridges(quantile_lines = TRUE) +
  xlim(40, 90) +
  labs(x = "\nLength (cm)", y = "Count\n") + 
  theme(legend.position = c(0.8, 0.7))

# Caption: Length frequency distribution by escape ring treatment.
ggsave(paste0("figures/size_freq_", YEAR, ".png"), dpi=300, height=3, width=6, units="in")

# The control pot captured 

