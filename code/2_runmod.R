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

# Priors for s10, s50, s90 parameterization
priors1 <- sprior %>% 
  mutate(p = round(p, 1)) %>% 
  filter(p %in% c(0.1, 0.5, 0.9) & Treatment != "Control") %>% 
  group_by(Treatment, p) %>% 
  dplyr::summarise(prior = mean(length)) %>% 
  mutate(log_prior = log(prior))

# s50 from theoretical selectivities
sprior %>% 
  mutate(p = round(p, 2)) %>% 
  filter(p %in% c(0.50) & Treatment != "Control") %>% 
  group_by(Treatment, p) %>% 
  dplyr::summarise(prior = mean(length)) %>% 
  ungroup() %>% 
  mutate(trt = 2.54 * c(3.5, 3.75, 4)) -> s50_thx

ggplot(s50_thx, aes(x = trt, y = prior)) +
  geom_point() +
  geom_smooth(method = "lm")

# Prior selection range (25-75%)
sprior %>% 
  mutate(p = round(p, 2)) %>%
  filter(p %in% c(0.25, 0.75) & Treatment != "Control") %>% 
  group_by(Treatment, p) %>% 
  dplyr::summarise(prior = mean(length)) %>% 
  ungroup() %>% 
  data.table::dcast(Treatment ~ p) %>% 
  mutate(SR = `0.75` - `0.25`) %>% 
  mutate(trt = 2.54 * c(3.5, 3.75, 4)) -> sr_thx

ggplot(sr_thx, aes(x = trt, y = SR)) +
  geom_point() +
  geom_smooth(method = "lm")

sc <- 8
sh <- 3
(mu <- sc / (sh - 1))
x=seq(from=0.01, to=20, by= 0.1)
data=dinvgamma(x, scale = sc, shape = sh)
plot(x, data)
# sprior %>% 
#   mutate(p = round(p, 1)) %>% 
#   filter(p == 0.5 & Treatment != "Control") %>% 
#   group_by(Treatment) %>% 
#   dplyr::summarise(s50 = mean(length)) %>% 
#   mutate(log_s90 = log(s50))

# Summarize counts

tagged <- counts %>% 
  group_by(set, Treatment, disposition) %>% 
  dplyr::summarise(n = n()) %>% 
  filter(disposition == "tagged")

# Remove outliers

bio %>% filter(length != min(bio$length)) %>% distinct(length) %>% arrange(length)
bio %>% filter(length == 30)
# Length bins
# Same length bin structure used by Feds)
# bio <- bio %>%
#   filter(!c(length < 40)) %>%
#   mutate(length2 = ifelse(length < 41, 41,
#                           ifelse(length > 99, 99, length)),
#          length_bin = cut(length2, breaks = seq(39.9, 99.9, 2),
#                           labels = paste(seq(41, 99, 2)))) %>%
#   select(-length2)

# Length bins used in BC study Haist et al 2004
# bio <- bio %>%
#   filter(!c(length < 50)) %>%
#   mutate(length2 = ifelse(length < 51, 51,
#                           ifelse(length > 95, 95, length)),
#          length_bin = cut(length2, breaks = seq(50.9, 94.9, 1),
#                           labels = paste(seq(51, 94, 1)))) %>%
#   select(-length2)
bio <- bio %>%
  filter(!c(length < 50)) %>%
  mutate(length2 = ifelse(length < 51, 51,
                          ifelse(length > 79, 79, length)),
         length_bin = cut(length2, breaks = seq(50.9, 79.9, 1),
                          labels = paste(seq(51, 79, 1)))) %>%
  select(-length2)

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

# ggplot(sum_df, aes(x = length_bin, y = p,
#                    group =  factor(effort_no), col = factor(effort_no))) +
#   geom_point() +
#   geom_line() +
#   facet_grid(effort_no ~ Treatment)

# remove_sets <- c(3,8,12,13,14,15,17)
# trt %>% group_by(effort_no) %>% summarise(sum(tot_n))
# trt <- trt %>% filter(!c(effort_no %in% remove_sets))

# trt <- mutate(trt, Effort_no = factor(effort_no))


# mod <- glm(p ~ length_bin, family = "binomial", weights = tot_n, data = trt)
# mod <- glmer(p ~ length_bin + (1 | effort_no), family = "binomial", weights = tot_n, data = trt)
# tidy(mod)
# summary(mod)
# pred_df <- data.frame(length_bin = seq(30, 100, 1)) %>% 
#   mutate(effort_no = 6)
# pred_df$pred <- predict(mod, pred_df, type = "response")

# Compile TMB code and fit model ----
setwd("~/great_escape/code")

# Escape 2 ----

# Third model attempt that links selectivity parameter estimates between
# treatments and has a global delta

df <- sum_df
# df <- com

# Number of pots per set that were sampled for each treatment
ntrt_df <- counts %>% 
  group_by(Treatment, set) %>% 
  dplyr::summarise(ntrt = length(which(disposition == "tagged"))) %>% 
  dcast(set ~ Treatment)

# Priors on 100 retention:
s100 <- filter(sprior, p == 1 & Treatment != "Control") %>% 
  group_by(Treatment) %>% 
  filter(length == min(length)) %>% 
  select(Treatment, length)
s0 <- filter(sprior, p < 0.01 & Treatment != "Control") %>%
  group_by(Treatment) %>% 
  filter(length == max(length)) %>% 
  select(Treatment, length)

len <- sort(unique(df$length_bin))

# TMB index for s0 and s100 prior length (TMB indexing starts with 0)
s0_index <- c(0, 3, 7)
s100_index <- c(20, 24, 27)

# Rescale lengths
sprior %>% 
  filter(length %in% unique(df$length_bin)) %>% 
  mutate(scaled = scale(length)[,1]) -> sprior

# Theoretical s50 and slope for each treatment
s50_vec <- vector(length = 3)
slp_vec <- vector(length = 3)
alpha_vec <- vector(length = 3)
beta_vec <- vector(length = 3)

for(i in 1:length(unique(df$Treatment))) {
  tmp <- sprior %>% filter(Treatment == unique(df$Treatment)[i])
  # fit <- glm(p ~ scaled, data = tmp, family = "quasibinomial")
  fit <- glm(p ~ length, data = tmp, family = "quasibinomial")
  s50_vec[i] <- - coef(fit)[1] / coef(fit)[2]
  slp_vec[i] <- coef(fit)[2] / 4 # FLAG - negative?
  alpha_vec[i] <- coef(fit)[1]
  beta_vec[i] <- coef(fit)[2]
}

plot(fitted(fit) ~ tmp$scaled)

# Adjust other quantities for scaling
sc_len <- scale(len)[,1]
mu <- mean(unique(df$length_bin))
sd <- sd(unique(df$length_bin))
fit_len <- seq(10, 100, 1)
sc_fit_len <- (mu - fit_len) / sd
# fit_len <- len
trt_len <- 2.54 * c(3.5, 3.75, 4) # diameter of escape ring treatments in cm
sc_trt_len <- (mu - trt_len) / sd
# mu - trt_len * sd # method to backtransform

# Regress on raw, not scaled
lm_s50 <- coef(lm(s50_vec ~ trt_len))
lm_slp <- coef(lm(slp_vec ~ trt_len))
plot(trt_len, s50_vec)
plot(trt_len, slp_vec)

data <- list(slx_type = 1, # model switch
             nset = length(unique(df$effort_no)),#17, # # number of sets
             nlen = length(unique(df$length_bin)), # number of length bins
             ntrt = length(unique(df$Treatment)), # number of treatments
             len = len, # vector of lengths for which there are data
             fit_len = seq(10, 100, 1), #fit_len, # vector of lengths for fitted values
             trt = trt_len, # vector of escape ring treatment diameters
             npot_ctl = ntrt_df$Control, # vector of number of control pots in each set
             npot_trt = as.matrix(select(ntrt_df, matches(paste(unique(df$Treatment), collapse = "|")))), # matrix of number of experimental pots in each set [j,k]
             ctl_dat = matrix(df$ctl_n[df$Treatment == "3.50 in"], ncol = length(unique(df$effort_no))), #1), # # control pots: matrix of number of fish caught by length bin (row) by set (col)
             trt_dat = array(df$exp_n, # experimental pots: matrix of number of fish caught by length bin (row) by set (col)
                             dim = c(length(unique(df$length_bin)),
                                     length(unique(df$effort_no)),
                                     length(unique(df$Treatment)))),
             theor_s50 = s50_vec, # theoretical s50s
             theor_slp = slp_vec, # theoretical slopes
             wt_s50 = 1,  # weight for selectivity penalties                            
             wt_slp = 1,
             sigma_s0 = 0.4,  # priors to constrain selectivity at 0 and 1
             sigma_s100 = 0.05,
             s0_index = s0_index,
             s100_index = s100_index) 

parameters <- list(dummy = 0,
                   # alpha = alpha_vec,
                   # beta = beta_vec,
                   s50 = s50_vec,
                   slp = slp_vec,
                   log_delta = log(1),
                   a1 = lm_s50[1],
                   log_b1 = log(lm_s50[2]),
                   # a2 = lm_slp[1],
                   # b2 = lm_slp[2],
                   nu = rep(0, length(unique(df$effort_no))))

# Map
map <- list(dummy = factor(NA),#
            # alpha = rep(factor(NA), 3),
            # beta = rep(factor(NA), 3),
            # log_delta = factor(NA),
            a1 = factor(NA),
            b1 = factor(NA),
            a2 = factor(NA),
            b2 = factor(NA)
            # nu = rep(factor(NA), length(unique(df$effort_no)))
)

compile("escape2.cpp")
dyn.load(dynlib("escape2"))

# Bounds
l_alpha <- rep(-45, length(unique(df$Treatment)));
u_alpha <- rep(-20, length(unique(df$Treatment)));

l_beta <- rep(0.4, length(unique(df$Treatment)));
u_beta <- rep(0.8, length(unique(df$Treatment)));

l_log_delta <- log(0.5);
u_log_delta <- log(1.5);

l_a1 <- 0.01;
u_a1 <- 5;
l_b1 <- 0.01;
u_b1 <- 10;
l_a2 <- 0.01;
u_a2 <- 5;
l_b2 <- -5;
u_b2 <- -0.01;

l_nu <- rep(-5, data$nset)
u_nu <- rep(5, data$nset)
# lowbnd <- c(l_alpha, l_beta, l_log_delta, l_a1, l_b1, l_a2, l_b2, l_nu) 
# uppbnd <- c(u_alpha, u_beta, u_log_delta, u_a1, u_b1, u_a2, u_b2, u_nu) 

# Model 1: Assume fixed delta
map <- list(dummy = factor(NA),
            log_delta = factor(NA))
map <- list(dummy = factor(NA))
lowbnd <- c(l_alpha, l_beta, l_a1, l_b1, l_a2, l_b2, l_nu) 
uppbnd <- c(u_alpha, u_beta, u_a1, u_b1, u_a2, u_b2, u_nu) 

# Model 2: Assume fixed delta and no random effects
# map <- list(dummy = factor(NA),#
#             log_delta = factor(NA),
#             nu = rep(factor(NA), length(unique(df$effort_no))))
# lowbnd <- c(l_alpha, l_beta, l_a1, l_b1, l_a2, l_b2) 
# uppbnd <- c(u_alpha, u_beta, u_a1, u_b1, u_a2, u_b2) 

model <- MakeADFun(data, parameters, map = map, 
                   DLL = "escape2", silent = TRUE,
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
     col = "black", xlab = "Length (cm)", ylab = "# Treatment / (# Control + # Treatment)", ylim = c(0,1))
plot(fit_len, model$report()$fit_phi[,1], type = "l", 
     col = "black", xlab = "Length (cm)", ylab = "# Treatment / (# Control + # Treatment)", ylim = c(0,1))

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

compile("escape.cpp")
dyn.load(dynlib("escape"))
compile("escape3.cpp")
dyn.load(dynlib("escape3"))

# TREATMENT <- "3.50 in"
# df = sum_df;
# treatment = TREATMENT; 
# model = MODEL;
# # Priors
# #mu_log_s90 = mu_log_s90;
# mu_log_s10 = log(60);
# sig_log_s90 = log(2);
# sig_log_s10 = log(2);
# # Initial values ("i_"): Starting values from Table N.4
# # Model 1.3 and 1.2 in Haist et al. 2004
# i_log_s50 = log(65);
# i_log_s90 = log(69);
# i_log_s10 = log(60);
# # Lower ("l_") and upper ("u_") bounds
# l_log_s50 = log(55);
# u_log_s50 = log(65);
# l_log_s90 = log(65);
# u_log_s90 = log(80);
# l_log_s10 = log(40);
# u_log_s10 = log(50);
# l_log_delta = log(0.5);
# u_log_delta = log(1.5) 

# Function to fit selectivity curves
fit_slx <- function(df = sum_df,
                    treatment = TREATMENT, 
                    model = MODEL,
                    # Priors
                    #mu_log_s90 = mu_log_s90,
                    mu_log_s10 = log(60),
                    sig_log_s90 = log(5),
                    sig_log_s10 = log(5),
                    # Initial values ("i_"): Starting values from Table N.4
                    # Model 1.3 and 1.2 in Haist et al. 2004
                    i_log_s50 = log(65),
                    i_log_s90 = log(69),
                    i_log_s10 = log(60),
                    # Lower ("l_") and upper ("u_") bounds
                    l_log_s50 = log(55),
                    u_log_s50 = log(65),
                    l_log_s90 = log(65),
                    u_log_s90 = log(75),
                    l_log_s10 = log(45),
                    u_log_s10 = log(50),
                    l_log_delta = log(0.5),
                    u_log_delta = log(1.5) 
                    ) {

  trt <- filter(df, Treatment == TREATMENT)  # Subset treatment
  # trt <- mutate(trt, effort_no = 1)
  
  # Default data and switches
  data <- list(slx_type = 1, # model switch
               nset = length(unique(trt$effort_no)), # number of sets
               nlen = length(unique(trt$length_bin)), # number of length bins
               len = sort(unique(trt$length_bin)), # vector of lengths for which there are data
               fitted_len = seq(10, 100, 1), # vector of lengths for fitted values
               npot_ctl = tagged %>% filter(Treatment == "Control") %>% pull(n), # vector of number of control pots in each set
               npot_exp = tagged %>% filter(Treatment == TREATMENT) %>% pull(n), # vector of number of experimental pots in each set
               ctl_dat = matrix(trt$ctl_n, ncol = length(unique(trt$effort_no))), # control pots: matrix of number of fish caught by length bin (row) by set (col)
               exp_dat = matrix(trt$exp_n, ncol = length(unique(trt$effort_no))), # experimental pots: matrix of number of fish caught by length bin (row) by set (col)
               # Priors
               mu_log_s90 = priors1 %>% filter(Treatment %in% TREATMENT & p == 0.9) %>% pull(log_prior), # prior on log s90 from theoretical selectivity curves
               mu_log_s10 = priors1 %>% filter(Treatment %in% TREATMENT & p == 0.1) %>% pull(log_prior), 
               sig_log_s90 = sig_log_s90,
               sig_log_s10 = sig_log_s10) 
  
  parameters <- list(dummy = 0,
                     log_s50 = i_log_s50,
                     log_s90 = i_log_s90,
                     log_s10 = i_log_s10,
                     log_delta = log(1),
                     nu = rep(0, length(unique(trt$effort_no))))
  
  # Model 0: Troubleshooting mode
  if (model == 0) 
  {
    map <- list(log_s50 = factor(NA),
                log_s90 = factor(NA),
                log_s10 = factor(NA),
                log_delta = factor(NA),
                nu = rep(factor(NA), length(unique(trt$effort_no))))
    lownd <- c(0)
    uppbnd <- c(0)
    
    obj <- MakeADFun(data, parameters, map = map, DLL = "escape3")
  }
  
  # Model 1: Fixed delta, symmetric logistic (2 param)     
  if (model == 1)
  {
    map <- list(dummy = factor(NA), 
                log_s10 = factor(NA),
                log_delta = factor(NA))#,
                # nu = rep(factor(NA), length(unique(trt$effort_no))))
    # Bound random effect between [-5, 5]
    lowbnd <- c(l_log_s50, l_log_s90, rep(-5, data$nset)) 
    uppbnd <- c(u_log_s50, u_log_s90, rep(5, data$nset)) 
    obj <- MakeADFun(data, parameters, map = map, 
              DLL = "escape3", silent = FALSE,
              hessian = TRUE, random = "nu")
  }
  
  # Model 2: Fixed delta, asymmetric logistic (3 param)
  if (model == 2)
  {
    data$slx_type <- 2
    map <- list(dummy = factor(NA),
                log_delta = factor(NA))#,
                # nu = rep(factor(NA), length(unique(trt$effort_no))))
    lowbnd <- c(l_log_s50, l_log_s90, l_log_s10, rep(-5, data$nset)) 
    uppbnd <- c(u_log_s50, u_log_s90, u_log_s10, rep(5, data$nset)) 
    obj <- MakeADFun(data, parameters, map = map, 
                     DLL = "escape3", silent = FALSE,
                     hessian = TRUE, random = "nu")
  }
  
  # Model 3: Estimate delta, symmetric logistic (3 param)
  if (model == 3)
  {
    map <- list(dummy = factor(NA), 
                log_s10 = factor(NA))#,
                # nu = rep(factor(NA), length(unique(trt$effort_no))))
    lowbnd <- c(l_log_s50, l_log_s90, l_log_delta, rep(-5, data$nset)) 
    uppbnd <- c(u_log_s50, u_log_s90, u_log_delta, rep(5, data$nset)) 
    obj <- MakeADFun(data, parameters, map = map, 
                     DLL = "escape3", silent = FALSE,
                     hessian = TRUE, random = "nu")
  }
  
  # Model 3: Estimate delta, asymmetric logistic (4 param)
  if (model == 4)
  {
    data$slx_type <- 2
    map <- list(dummy = factor(NA))#,
                # nu = rep(factor(NA), length(unique(trt$effort_no))))
    lowbnd <- c(l_log_s50, l_log_s90, l_log_s10, l_log_delta, rep(-5, data$nset)) 
    uppbnd <- c(u_log_s50, u_log_s90, u_log_s10, u_log_delta, rep(5, data$nset)) 
    obj <- MakeADFun(data, parameters, map = map, 
                     DLL = "escape3", silent = FALSE,
                     hessian = TRUE, random = "nu")
  }
  
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                lower = lowbnd, upper = uppbnd)
  rep <- sdreport(obj)
  out <- list("data" = data, "obj" = obj, "opt" = opt, "rep" = rep)
  return(out)
}

# User inputs
TREATMENT <- "4.00 in"
TREATMENT <- "3.75 in"
TREATMENT <- "3.50 in"

# Model options
# 0: Troubleshooting mode
# 1: Fixed delta, symmetric logistic
# 2: Fixed delta, asymmetric logistic
# 3: Estimate delta, symmetric logistic
# 4: Estimate delta, asymmetric logistic
MODEL <- 2

out <- fit_slx()

out[[3]]$convergence
out[[4]]
summary(out[[4]])
summary(out[[4]], "all", pvalue = TRUE)

plot(out[[1]]$fitted_len, out[[2]]$report()$fit_slx, type = "l", col = "black",
     xlab = "Length (cm)", ylab = "Selectivity")
if(out[[1]]$slx_type == 1) 
{
  abline(h = c(0.5, 0.9), 
         v = c(out[[2]]$report()$s50, out[[2]]$report()$s90),
         col = "grey", lty = 2)
} else {
  abline(h = c(0.1, 0.5, 0.9), 
         v = c(out[[2]]$report()$s10, out[[2]]$report()$s50, out[[2]]$report()$s90),
         col = "grey", lty = 2)
}


plot(out[[1]]$len, out[[2]]$report()$fit_phi, type = "l", col = "black",
     xlab = "Length (cm)", ylab = "# Treatment / (# Control + # Treatment)", ylim = c(0,1))
phi <- as.data.frame(out[[2]]$report()$phi)
names(phi) <- paste("set", 1:17, sep = "_")
for(i in 1:17) {
  tmp <- phi[,i]
  points(out[[1]]$len, tmp, add = TRUE, type = "l", col = "grey")
}
phi_obs <- rowSums(out[[1]]$exp_dat) / (rowSums(out[[1]]$exp_dat) + rowSums(out[[1]]$ctl_dat))

points(out[[1]]$len, phi_obs, add = TRUE)
points(out[[1]]$len, out[[2]]$report()$fit_phi, type = "l", add = TRUE)

sel <- array(data = NA,
             dim = c(length(out[[1]]$fitted_len), # nrow = length bins
                     length(unique(sum_df$Treatment)) + 1, # ncol = Treatments + row id
                     4)) # array dim = different model types

phi <- array(data = NA,
             dim = c(length(out[[1]]$fitted_len), # nrow = length bins
                     length(unique(sum_df$Treatment)) + 1, # ncol = Treatments + row id
                     4)) # array dim = different model types
pars <- list()
convergence <- list()

# Run escape2
compile("escape2.cpp")
dyn.load(dynlib("escape2"))

TREATMENT <- "3.50 in"
trt <- com %>% 
  filter(Treatment %in% TREATMENT) %>% 
  mutate(p = combined_p)
# trt <- filter(sum_df, Treatment == TREATMENT)  # Subset treatment
trt <- mutate(trt, effort_no = 1)

data <- list(slx_type = 1, # model switch
             nset = length(unique(trt$effort_no)), # number of sets
             nlen = length(unique(trt$length_bin)), # number of length bins
             len = sort(unique(trt$length_bin)), # vector of lengths for which there are data
             fitted_len = seq(10, 100, 1), # vector of lengths for fitted values
             npot_ctl = rep(4, length(unique(trt$effort_no))), # vector of number of control pots in each set
             npot_exp = rep(4, length(unique(trt$effort_no))), # vector of number of experimental pots in each set
             ctl_dat = matrix(trt$ctl_n, ncol = length(unique(trt$effort_no))), # control pots: matrix of number of fish caught by length bin (row) by set (col)
             exp_dat = matrix(trt$exp_n, ncol = length(unique(trt$effort_no)))) # experimental pots: matrix of number of fish caught by length bin (row) by set (col)

# data.frame(x = seq(10, 100, 2),
#            scaled = scale(seq(10, 100, 2)))
a <- 13; b <- -0.22
lengths <- seq(10, 100, 2)
# lengths <- scale(seq(10, 100, 2))
slx <- vector(length = length(lengths))
for(i in 1:length(lengths)) {
  # slx[i] <- 1 / (1 + exp(a + b * lengths[i]))
  # slx[i] <- exp(a + b * lengths[i]) / (1 + exp(a + b * lengths[i]))
  slx[i] <- 1 / (1 + exp(-(a + b * lengths[i])))
}
plot(lengths, slx)

parameters <- list(dummy = 0,
                   alpha = a,
                   beta = b,
                   log_delta = log(1),
                   nu = rep(0, length(unique(trt$effort_no))))

map <- list(dummy = factor(NA),
            # alpha = factor(NA),
            log_delta = factor(NA))#,
            # nu = rep(factor(NA), length(unique(trt$effort_no))))

model <- MakeADFun(data, parameters, map = map, 
                   DLL = "escape2", silent = FALSE,
                   hessian = TRUE, random = "nu") #
# checking for minimization
xx <- model$fn(model$env$last.par)
print(model$report())

fit <- nlminb(model$par, model$fn, model$gr, 
              # control=list(eval.max=1000000,iter.max=100000),
              control=list(rel.tol=1e-12,
                           eval.max=100000,iter.max=10000))#,
              # lower=lowbnd,upper=uppbnd)
best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(rep)
plot(data$fitted_len, model$report()$fit_slx, type = "l", col = "black", ylim = c(0,1))
slx <- as.data.frame(model$report()$slx)
names(slx) <- paste("set", 1:17, sep = "_")
for(i in 1:17) {
  tmp <- slx[,i]
  points(data$len, tmp, add = TRUE, type = "l", col = "grey")
}

plot(data$fitted_len, model$report()$fit_phi, ylim = c(0,1), type = "l", col = "black")
phi <- as.data.frame(model$report()$phi)
names(phi) <- paste("set", 1:17, sep = "_")
for(i in 1:17) {
  tmp <- phi[,i]
  points(data$len, tmp, add = TRUE, type = "l", col = "grey")
}

#-------

#------
(s50 <- model$report()$s50)
(s10 <- model$report()$s10)
(s90 <- model$report()$s90)
(fit_slx <- model$report()$fit_slx)
# r <- exp(best[2])# model$report()$slxr
(phi <- model$report()$phi)

# slx <- vector(length = data$nlen)
lengths <- seq(10, 100, 2)
slx <- vector(length = length(lengths))
# for(i in 1:data$nlen) {
#   len <- data$len[i]
for(i in 1:length(lengths)) {
  len <- lengths[i]
  
  slx[i] <- 
  # if(data$model == 1) {
  #   if(len <= s50) {
  #     delta_slx <- s50 - (2 * s50 - s90)
  #   } else {
  #     delta_slx <- s90 - s50
  #   }
  # } else {
  #   if(len <= s50) {
  #     delta_slx <- s50 - s10
  #   } else {
  #     delta_slx <- s90 - s50
  #   }
  # }
  # slx[i] <- 1 / (1 + exp(-2 * log(3) * ((len - s50) / delta_slx)))

plot(lengths, slx, col = "red")
points(data$fitted_len, fit_slx, type = "l", col = "black")



# -----
cat(model$report()$pfit,"\n")
res <- data.frame(p = exp(model$report()$log_s50),
           p_est = rep(as.list(rep, what = "Estimate")$`p`, length(model$report()$p)),
           p_std = rep(as.list(rep, what = "Std")$`p`, length(model$report()$p)))

ggplot(data = res, aes(x = p)) +
  geom_histogram(bins = 20, alpha = 0.5, aes(x = p, y = ..density..), position = "identity") +
  geom_vline(xintercept = res$p_est) +
  geom_vline(xintercept = res$p_est + (1.96 * res$p_std), lty = 3) +
  geom_vline(xintercept = res$p_est - (1.96 * res$p_std), lty = 3) +
  ggtitle(label = "Model 1") +
  theme(plot.title = element_text(hjust = 0)) -> mod1

ggsave("mod1_fit.png", plot = mod1, dpi = 300, height = 2.5, width = 5, units = "in")

# 4b. Model 2 ----

data$model <- 2
map <- list(dummy = factor(NA))
model <- MakeADFun(data, parameters, map=map, 
                   DLL="HW2_jysullivan",silent=T,
                   hessian=T)

# checking for minimization
xx <- model$fn(model$env$last.par)
print(model$report())
# Bounds on the parameters - needs to be equal to the number of ESTIMATED
# parameters not equal to the length of the total number of pars in the
# parameters vectors (dummy disappears) use:
length(model$par) # how many parameters for upper and lower bounds
lowbnd= c(0, # p
          0) # tau

uppbnd= c(1, # p
          Inf) # tau

fit <- nlminb(model$par, model$fn, model$gr, 
              control=list(rel.tol=1e-12,
                           eval.max=100000,iter.max=100000),
              lower=lowbnd,upper=uppbnd)

best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(rep)
cat(model$report()$pfit,"\n")
cat(model$report()$log_sigma_y,"\n")
model$report()$p
as.list(rep, what = "Std")$`p`

res <- data.frame(p = model$report()$p,
                  p_est = rep(as.list(rep, what = "Estimate")$`p`, length(model$report()$p)),
                  p_std = rep(as.list(rep, what = "Std")$`p`, length(model$report()$p)))

ggplot(data = res, aes(x = p)) +
  geom_histogram(bins = 20, alpha = 0.5, aes(x = p, y = ..density..), position = "identity") +
  geom_vline(xintercept = res$p_est) +
  geom_vline(xintercept = res$p_est + (1.96 * res$p_std), lty = 3) +
  geom_vline(xintercept = res$p_est - (1.96 * res$p_std), lty = 3) +
  ggtitle(label = "Model 2") +
  theme(plot.title = element_text(hjust = 0)) -> mod2

ggsave("mod2_fit.png", plot = mod2, dpi = 300, height = 2.5, width = 5, units = "in")

# 3. Simulation ----

set.seed(222)
sim_sexratio <- function(reps = 1000, # number of simulations
                         nyr = 25, # number of years 
                         nsamps = 100, # number of samples collected each yr)
                         mod = 1) { # which model to use
  
  # True p (expected value for p when p ~ beta(2,1))
  p_true <- mean(rbeta(reps*nyr, 2, 1))
  
  # Different p for each year and simulation 
  dat_ls <- list()
  MM <- vector()
  FF <- vector()
  
  for(i in 1:reps){
    # simulated p's, where p~beta(2,1)
    yr_p <- rbeta(nyr, 2, 1)
    
    for (j in 1:nyr) {
      # simulate data using the binomial distribution and our known p
      MM[j] <- rbinom(n = 1, size = nsamps, prob = yr_p[j])
      FF[j] <- 100 - MM[j]
    }
    dat_ls[[i]] <- list(MM = MM, FF = FF)
  }
  
  results <- list()
  bool <- vector()
  year <- vector()
  
  # mod = 2
  for(i in 1:reps) {
    
    data <- dat_ls[[i]]
    
    if(mod == 1) { 
      data$model <- 1
    } else {
      data$model <- 2
    } 
    
    model <- MakeADFun(data, parameters, map=map, 
                       DLL="HW2_jysullivan",silent=T,hessian=T)
    
    fit <- nlminb(model$par, model$fn, model$gr, 
                  control=list(rel.tol=1e-12,
                               eval.max=100000,iter.max=10000),
                  lower=lowbnd,upper=uppbnd)
    
    best <- model$env$last.par.best
    rep <- sdreport(model)
    
    p_est <- as.list(rep, what = "Est")$`p`
    p_se <- as.list(rep, what = "Std")$`p`
    p_lower95 <- p_est - 1.96 * p_se
    p_upper95 <- p_est + 1.96 * p_se
    
    for(j in 1:nyr) {
      bool[j] <- ifelse(p_true <= p_upper95 & p_true >= p_lower95, 1, 0)
      year[j] <- j
    }
    
    results[[i]] <- data.frame(p_est = rep(p_est, nyr),
                               p_se = rep(p_se, nyr),
                               p_lower95 = rep(p_lower95, nyr),
                               p_upper95 = rep(p_upper95, nyr),
                               p_true = rep(p_true, nyr), bool = bool,
                               year = year, sim = i)
    
  }
  return(results)
}

# Sim Model 1 ----
lowbnd= c(0) # p
uppbnd= c(1) # p
map <- list(dummy = factor(NA), tau = factor(NA))
results <- sim_sexratio(reps = 1000, # number of simulations
                        nyr = 25, # number of years 
                        nsamps = 100, # number of samples collected each yr)
                        mod = 1)
res <- do.call("rbind", results)

res %>% 
  group_by(bool) %>% 
  summarize(Count = n()) %>% 
  mutate(bool = ifelse(bool==0, "No", "Yes"),
         Percent = formatC(round(Count/length(res$p_est)*100, 1), small.interval = 1)) %>% 
  select(`95% CI contains true p?` = bool, Count, Percent) %>% 
  kable()

res %>% 
  ggplot(aes(x = p_est, fill = factor(bool), colour = factor(bool))) + 
  geom_histogram(bins = 100, alpha = 0.5, 
                 aes(y = ..density.., fill = factor(bool), colour = factor(bool)), 
                 position = "identity") +
  geom_vline(xintercept = res$p_true, size = 1, lty = 2) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  scale_colour_grey(start = 0.3, end = 0.8) +
  theme(legend.position = "none") + 
  ggtitle(label = "Model 1 simulation") +
  labs(x = "Estimated p", y = "Density") +
  theme(plot.title = element_text(hjust = 0)) -> psims

ggsave("psims_model1.png", plot = psims, dpi = 300, height = 2.5, width = 5, units = "in")

# Sim Model 2 ----
map <- list(dummy = factor(NA))
lowbnd= c(0, # p
          0) # tau

uppbnd= c(1, # p
          Inf) # tau
results <- sim_sexratio(reps = 1000, # number of simulations
                        nyr = 25, # number of years 
                        nsamps = 100, # number of samples collected each yr)
                        mod = 2)
res <- do.call("rbind", results)

res %>% 
  group_by(bool) %>% 
  summarize(Count = n()) %>% 
  mutate(bool = ifelse(bool==0, "No", "Yes"),
         Percent = formatC(round(Count/length(res$p_est)*100, 1), small.interval = 1)) %>% 
  select(`95% CI contains true p?` = bool, Count, Percent) %>% 
  kable()

res %>% 
  ggplot(aes(x = p_est, fill = factor(bool), colour = factor(bool))) + 
  geom_histogram(bins = 70, alpha = 0.5, 
                 aes(y = ..density.., fill = factor(bool), colour = factor(bool)), 
                 position = "identity") +
  geom_vline(xintercept = res$p_true, size = 1, lty = 2) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  scale_colour_grey(start = 0.3, end = 0.8) +
  theme(legend.position = "none") + 
  ggtitle(label = "Model 2 simulation") +
  labs(x = "Estimated p", y = "Density") +
  theme(plot.title = element_text(hjust = 0))-> psims

ggsave("psims_model2.png", plot = psims, dpi = 300, height = 2.5, width = 5, units = "in")

dat
ggplot(dat, aes(x = Year, y = pM)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  ylab("Proportion male") -> ts

ggsave("ts.png", plot = ts, dpi = 300, height = 2.5, width = 5, units = "in")

acf(dat$pM)
