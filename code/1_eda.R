# Exploratory data analysis
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-10-01

# Set up ----

YEAR <- 2019 # study year(s)
source("code/helper.r")
bio <- read_csv(paste0("data/pot_bio_", YEAR, ".csv")) %>% 
  filter(!is.na(length)) %>% 
  mutate(Treatment = derivedFactor("Control" = Treatment == "99",
                                   "3.50 in" = Treatment == "02",
                                   "3.75 in" = Treatment == "01",
                                   "4.00 in" = Treatment == "00",
                                   .default = NA,
                                   .ordered = TRUE))

# Email sent to A. Baldwin 2019-10-01 about the NAs. FIX ME
bio %>% filter(is.na(Treatment))

# Sample sizes
bio %>% 
  group_by(Treatment) %>% 
  filter(!is.na(Treatment)) %>% 
  dplyr::summarise(n = n(),
            mean = mean(length),
            median = median(length)) %>% 
  kable()

# Number of pots
n_distinct(bio$effort_no)

# Length comp ggridge plots
bio %>% 
  filter(!is.na(Treatment) & between(length, 40, 90)) %>% 
  droplevels() %>% 
  ggplot(aes(x = length, y = Treatment, group = Treatment, fill = Treatment)) +
  geom_density_ridges(aes(point_fill = Treatment, point_color = Treatment),
                      alpha = 0.3) +
  scale_fill_grey() +
  # stat_density_ridges(quantile_lines = TRUE) +
  xlim(40, 90) +
  labs(x = "\nLength (cm)", y = NULL) + 
  ggtitle("Length compositions by escape ring treatment") +
  # scale_y_reverse() +
  theme(legend.position = "none")


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

# Caption: Size frequency distribution by escape ring treatment.
ggsave(paste0("figures/size_freq_", YEAR, ".png"), dpi=300, height=3, width=6, units="in")


# The control pot captured 

remotes::install_github("mebrooks/selfisher/selfisher", build_vignette = TRUE)
install.packages("glmmTMB")
install.packages("bbmle")
library(selfisher)
library(TMB)
library(plyr)
library(ggplot2); theme_set(theme_bw())
data("ccmhsdat")
head(ccmhsdat)

## ----aggregate-----------------------------------------------------------

sumdat=ddply(ccmhsdat, ~length+type, summarize, prop=sum(test)/sum(total), total=sum(total))

ccmhsdat %>% group_by(type, haul) %>% dplyr::summarise(sum(prop))
## ----fits----------------------------------------------------------------
mod_both=selfisher(prop~length*type, total=total, ccmhsdat, haul=haul)

## ----pred----------------------------------------------------------------
newdata=expand.grid(length=unique(ccmhsdat$length),
                    total=1,
                    haul=1,
                    type=c("baseline", "stimulation"))

newdata$prop=predict(mod_both, newdata=newdata, type="response")

## ----ci------------------------------------------------------------------
bs=bootSel(mod_both, nsim=100, parallel = "multicore", ncpus = 4, FUN=function(mod){predict(mod, newdata=newdata, type="response")})

quants=apply(bs$t, 2, quantile, c(0.025, 0.5, 0.975))
newdata[,c("lo", "mid", "hi")]=t(quants)

## ----plot----------------------------------------------------------------
ggplot(sumdat, aes(length, prop, colour=type))+geom_point(aes(size=total), alpha=0.5)+
  geom_line(data=newdata)+
  geom_ribbon(data=newdata, aes(ymin=lo, ymax=hi, fill=type), alpha=0.2)
