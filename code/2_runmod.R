# Simple logistic regression
# Jane Sullivan
# jane.sullivan1@alaska.gov
# Last updated 2019-10-07

source("code/helper.r")
library(TMB)
library(cowplot)
library(nimble) # inverse gamma prior
library(broom) # tidy() clean up parameter estimates

# Data ----
YEAR <- 2019
bio <- read_csv(paste0("data/pot_bio_", YEAR, ".csv")) %>% 
  filter(!is.na(length)) %>% 
  mutate(Treatment = derivedFactor("Control" = Treatment == "99",
                                   "3.50 in" = Treatment == "02",
                                   "3.75 in" = Treatment == "01",
                                   "4.00 in" = Treatment == "00",
                                   .default = NA,
                                   .ordered = TRUE))
# Remove data that don't have treatment data associated with them
bio <- bio %>% filter(!is.na(Treatment))

# Total sablefish counts by depth and disposition (some pots in a set were
# dumped due to processing time)
counts <- read_csv(paste0("data/total_counts_", YEAR, ".csv")) %>% 
  mutate(Treatment = derivedFactor("Control" = treatment == "Blue",
                            "3.50 in" = treatment == "Purple",
                            "3.75 in" = treatment == "Green",
                            "4.00 in" = treatment == "Yellow",
                            .default = NA,
                            .ordered = TRUE))

# Selectivity priors based on theoretical selectivity curves
sprior <- read_csv(paste0("output/theoretical_selectivity_", YEAR, ".csv"))

# Summarize counts
tagged <- counts %>% 
  group_by(set, Treatment, disposition) %>% 
  dplyr::summarise(n = n()) %>% 
  filter(disposition == "tagged")


# Length bins

# Bin structure used in Haist et al 2004. Examine lengths for which we believe
# the control pot is 100% selected
bio <- bio %>% filter(!c(length < 50)) %>%
mutate(length2 = ifelse(length < 51, 51, ifelse(length > 79, 79, length)),
length_bin = cut(length2, breaks = seq(50.9, 79.9, 1), labels = paste(seq(51,
79, 1)))) %>% select(-length2)

# Reorganize data and get number of fish caught in experimental pots / control
# pots for each length bin
sum_df <- bio %>% 
  group_by(Treatment, effort_no, length_bin, .drop=FALSE) %>% # 
  dplyr::summarise(n = n()) %>% 
  ungroup()

sum_df <- sum_df %>%
  filter(Treatment != "Control") %>% 
  dplyr::rename(exp_n = n) %>% 
  left_join(sum_df %>% 
              filter(Treatment == "Control") %>% 
              select(length_bin, effort_no, ctl_n = n)) %>% #
  mutate(tot_n = exp_n + ctl_n,
         # proportion retained in experimental out of control
         p = exp_n / tot_n,
         length_bin = as.numeric(as.character(length_bin)))

# show P combined
com <- bio %>% 
  group_by(Treatment, length_bin, .drop=FALSE) %>% # 
  dplyr::summarise(n = n()) %>% 
  ungroup()

com <- com %>%
  filter(Treatment != "Control") %>% 
  dplyr::rename(exp_n = n) %>% 
  left_join(com %>% 
              filter(Treatment == "Control") %>% 
              select(length_bin, ctl_n = n)) %>% #
  mutate(tot_n = exp_n + ctl_n,
         # proportion retained in experimental out of control
         combined_p = exp_n / tot_n,
         length_bin = as.numeric(as.character(length_bin)))

ggplot(sum_df %>% filter, aes(x = length_bin, y = p,
                   group =  Treatment, col = Treatment)) +
  geom_point() +
  facet_wrap(~Treatment) +
  ylim(c(0,1)) +
  geom_smooth()

ggplot(com, aes(x = length_bin, y = combined_p,
                   group =  Treatment, col = Treatment)) +
  geom_point() +
  ylim(c(0,1)) +
  geom_smooth()

# Model ----

setwd("~/great_escape/code")
df <- sum_df

# Number of pots per set that were sampled for each treatment
ntrt_df <- counts %>% 
  group_by(Treatment, set) %>% 
  dplyr::summarise(ntrt = length(which(disposition == "tagged"))) %>% 
  dcast(set ~ Treatment)

# Priors on 0 and 100% retention:
filter(sprior, p == 1 & Treatment != "Control") %>% 
  group_by(Treatment) %>% 
  filter(length == min(length)) %>% 
  select(Treatment, length)

filter(sprior, p < 0.01 & Treatment != "Control") %>%
  group_by(Treatment) %>% 
  filter(length == max(length)) %>% 
  select(Treatment, length)

len <- sort(unique(df$length_bin))
fit_len <- seq(10, 100, 1)

# TMB index for s0 and s100 prior length (TMB indexing starts with 0)
s0_index <- c(0, 3, 7)
s100_index <- c(20, 24, 27)

# Theoretical s50 and slope for each treatment
s50_vec <- vector(length = 3)
slp_vec <- vector(length = 3)

for(i in 1:length(unique(df$Treatment))) {
  tmp <- sprior %>% filter(Treatment == unique(df$Treatment)[i])
  fit <- glm(p ~ length, data = tmp, family = "quasibinomial")
  s50_vec[i] <- - coef(fit)[1] / coef(fit)[2]
  slp_vec[i] <- coef(fit)[2] / 4 
}

data <- list(nset = length(unique(df$effort_no)),#17, # # number of sets
             nlen = length(unique(df$length_bin)), # number of length bins
             ntrt = length(unique(df$Treatment)), # number of treatments
             len = len, # vector of lengths for which there are data
             fit_len = seq(10, 100, 1), #fit_len, # vector of lengths for fitted values
             npot_ctl = ntrt_df$Control, # vector of number of control pots in each set
             npot_trt = as.matrix(select(ntrt_df, matches(paste(unique(df$Treatment), collapse = "|")))), # matrix of number of experimental pots in each set [j,k]
             ctl_dat = matrix(df$ctl_n[df$Treatment == "3.50 in"], ncol = length(unique(df$effort_no))), #1), # # control pots: matrix of number of fish caught by length bin (row) by set (col)
             trt_dat = array(df$exp_n, # experimental pots: matrix of number of fish caught by length bin (row) by set (col)
                             dim = c(length(unique(df$length_bin)),
                                     length(unique(df$effort_no)),
                                     length(unique(df$Treatment)))),
             sigma_s0 = 0.35,  # priors to constrain selectivity at 0 and 1
             sigma_s100 = 0.05,
             s0_index = s0_index, # Index of length where theoretical curves were 0 and 1 for each treatment
             s100_index = s100_index) 

parameters <- list(dummy = 0,
                   s50 = s50_vec,
                   slp = slp_vec,
                   log_delta = log(1),
                   nu = rep(0, length(unique(df$effort_no))))

compile("escape.cpp")
dyn.load(dynlib("escape"))

# Bounds
l_log_delta <- log(0.5);
u_log_delta <- log(1.5);

l_nu <- rep(-5, data$nset)
u_nu <- rep(5, data$nset)
# lowbnd <- c(l_alpha, l_beta, l_log_delta, l_a1, l_b1, l_a2, l_b2, l_nu) 
# uppbnd <- c(u_alpha, u_beta, u_log_delta, u_a1, u_b1, u_a2, u_b2, u_nu) 

# Model 1: Assume fixed delta
map <- list(dummy = factor(NA),
            log_delta = factor(NA))
lowbnd <- c(l_alpha, l_beta, l_a1, l_b1, l_a2, l_b2, l_nu) 
uppbnd <- c(u_alpha, u_beta, u_a1, u_b1, u_a2, u_b2, u_nu) 

# Model 2: Estimate delta
map <- list(dummy = factor(NA))
# lowbnd <- c(l_alpha, l_beta, l_a1, l_b1, l_a2, l_b2) 
# uppbnd <- c(u_alpha, u_beta, u_a1, u_b1, u_a2, u_b2) 

model <- MakeADFun(data, parameters, map = map, 
                   DLL = "escape", silent = TRUE,
                   hessian = TRUE, random = "nu") #

# checking for minimization
xx <- model$fn(model$env$last.par)
print(model$report())

fit <- nlminb(model$par, model$fn, model$gr)#,  
# lower=lowbnd,upper=uppbnd)

best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(rep)

phi <- as.data.frame(model$report()$fit_phi)
names(phi) <- paste0(unique(df$Treatment))
phi <- phi %>% mutate(length_bin = len)#sort(unique(df$length_bin)))
phi <- melt(data = phi, id.vars = "length_bin", variable.name = "Treatment", value.name = "phi")

com %>% 
  left_join(phi, by = c("Treatment", "length_bin")) %>% 
  mutate(resid = combined_p - phi) -> tmp
  
ggplot(tmp) +
  geom_point(aes(x = length_bin, y = combined_p,
                      group =  Treatment, col = Treatment)) +
  # geom_hline(yintercept = 0.5) +
  geom_line(aes(x = length_bin, y = phi, group = Treatment, col = Treatment)) #+

ggplot(tmp, aes(x = length_bin, y = resid)) +
  geom_hline(yintercept = 0, col = "grey", size = 1) +
  geom_segment(aes(x = length_bin, xend = length_bin, y = 0, yend = resid),
               size = 0.2, col = "grey") +
  geom_point() +
  facet_wrap(~Treatment)

slx <- as.data.frame(model$report()$full_slx)
names(slx) <- paste0(unique(df$Treatment))
slx <- slx %>% mutate(length_bin = fit_len)# sort(unique(df$length_bin)))
slx <- melt(data = slx, id.vars = "length_bin", variable.name = "Treatment", value.name = "slx")

ggplot() +
  ylim(c(0,1)) +
  geom_line(data = slx, aes(x = length_bin, y = slx, group = Treatment, col = Treatment))

plot(sort(unique(df$length_bin)), model$report()$fit_phi[,1], type = "l", 
     col = "black", xlab = "Length (cm)", ylab = "# Treatment / (# Control + # Treatment)", 
     ylim = c(0,0.6))
phi <- as.data.frame(model$report()$phi[,,1])
names(phi) <- paste("set", 1:17, sep = "_")
for(i in 1:17) {
  tmp <- phi[,i]
  points(sort(unique(df$length_bin)), tmp, add = TRUE, type = "l", col = "grey")
}
model$report()$penl_s50
model$report()$penl_slp
model$report()$prior_s0
model$report()$prior_s100