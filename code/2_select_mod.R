# Escape ring selectivity model
# Jane Sullivan
# jane.sullivan1@alaska.gov
# Last updated 2019-11-26

source("code/helper.r")
library(TMB)
# library(nimble) # inverse gamma prior
# library(tmbstan) # run mcmc
# library(shinystan) # awesome diagnostic tool
set.seed(907)

# Data ----
YEAR <- 2019
bio <- read_csv(paste0("data/bio_cleaned_", YEAR, ".csv")) %>% 
  mutate(Treatment = derivedFactor("Control" = Treatment == "Control",
                                   "8.9 cm" = Treatment == "3.50 in",
                                   "9.5 cm" = Treatment == "3.75 in",
                                   "10.2 cm" = Treatment == "4.00 in",
                                   .default = NA,
                                   .ordered = TRUE))

# Total sablefish counts by depth and disposition (some pots in a set were
# dumped due to processing time)
counts <- read_csv(paste0("data/total_counts_", YEAR, ".csv")) %>% 
  mutate(Treatment = derivedFactor("Control" = treatment == "Blue",
                            "8.9 cm" = treatment == "Purple",
                            "9.5 cm" = treatment == "Green",
                            "10.2 cm" = treatment == "Yellow",
                            .default = NA,
                            .ordered = TRUE))

# Selectivity priors based on theoretical selectivity curves
sprior <- read_csv(paste0("output/theoretical_selectivity_", YEAR, ".csv")) %>%
  mutate(Treatment = derivedFactor("Control" = Treatment == "Control",
                                   "8.9 cm" = Treatment == "8.9 cm",
                                   "9.5 cm" = Treatment == "9.5 cm",
                                   "10.2 cm" = Treatment == "10.2 cm",
                                   .default = NA,
                                   .ordered = TRUE))

# Summarize counts
tagged <- counts %>% 
  group_by(set, Treatment, disposition) %>% 
  dplyr::summarise(n = n()) %>% 
  filter(disposition == "tagged")

# Length bins

# Bin structure used in Haist et al 2004. Examine lengths for which we believe
# the control pot is 100% selected. Bins are defined such that anything less than 
bio <- bio %>% 
  filter(!c(length < 50)) %>%
  mutate(length2 = ifelse(length < 51, 51, ifelse(length > 79, 79, length)),
         length_bin = cut(length2, breaks = seq(50.9, 79.9, 1), 
                          labels = paste(seq(51, 79, 1)))) %>% 
  select(-length2)

# Reorganize data and get number of fish caught in experimental pots / control
# pots for each length bin
df <- bio %>% 
  group_by(Treatment, effort_no, length_bin, .drop=FALSE) %>% # 
  dplyr::summarise(n = n()) %>% 
  ungroup()

df <- df %>%
  filter(Treatment != "Control") %>% 
  dplyr::rename(exp_n = n) %>% 
  left_join(df %>% 
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

ggplot(df %>% filter, aes(x = length_bin, y = p,
                   group =  Treatment, col = Treatment)) +
  geom_point() +
  facet_wrap(~Treatment) +
  ylim(c(0, 1)) +
  geom_smooth()

ggplot(com, aes(x = length_bin, y = combined_p,
                   group =  Treatment, col = Treatment)) +
  geom_point() +
  ylim(c(0,1)) +
  geom_smooth()

# Priors ---- 

# Priors on 0 and 100% retention. Include bin definitions to make indexing a
# little easier
prior_100 <- filter(sprior, p == 1 & Treatment != "Control") %>% 
  group_by(Treatment) %>% 
  filter(length == min(length)) %>% 
  select(Treatment, length) %>%
  mutate(length2 = ifelse(length < 51, 51, ifelse(length > 79, 79, length)),
         length_bin = cut(length2, breaks = seq(50.9, 79.9, 1), 
                          labels = paste(seq(51, 79, 1)))) %>% 
  select(-length2)
prior_100 # lengths at 100% retention in theoretical curves

prior_0 <- filter(sprior, p < 0.01 & Treatment != "Control") %>% 
  group_by(Treatment) %>% 
  filter(length == max(length)) %>% 
  select(Treatment, length) %>% 
  filter(!c(length < 50)) %>%
  mutate(length2 = ifelse(length < 51, 51, ifelse(length > 79, 79, length)),
         length_bin = cut(length2, breaks = seq(50.9, 79.9, 1), 
                          labels = paste(seq(51, 79, 1)))) %>% 
  select(-length2)
prior_0 # lengths at 0% retention in theoretical curves

# Index for length bins in TMB
tmb_index <- data.frame(length_bin = sort(unique(df$length_bin)),
                        index = 0:(length(sort(unique(df$length_bin)))-1))

len <- sort(unique(df$length_bin))
fit_len <- seq(10, 100, 1)

# TMB index for s0 and s100 prior length (TMB indexing starts with 0)
prior_0 <-  prior_0 %>% 
  mutate(length_bin = as.numeric(as.character(length_bin))) %>% 
  left_join(tmb_index)
s0_index <- prior_0 %>% pull(index)

prior_100 <- prior_100 %>% 
  mutate(length_bin = as.numeric(as.character(length_bin))) %>% 
  left_join(tmb_index) 
s100_index <- prior_100 %>% pull(index)

# Get starting values from theoretical selectivity curves
s50_vec <- vector(length = 3)
slp_vec <- vector(length = 3)

for(i in 1:length(unique(df$Treatment))) {
  tmp <- sprior %>% filter(Treatment %in% unique(df$Treatment)[i])
  fit <- glm(p ~ length, data = tmp, family = "quasibinomial")
  s50_vec[i] <- - coef(fit)[1] / coef(fit)[2]
  slp_vec[i] <- coef(fit)[2] / 4 
}

# Beta prior for s0 and s100
prior <- function(a, b, label){
  x <- seq(0.0001, .9999, 0.0001)
  val <- dbeta(x, a, b)
  return(data.frame('x' = x,
                    'y' = val,
                    'label'= rep(label, length(x))))
}

s0_alpha <- 1 # (0,1] to keep prior on border with 0; smaller values = tighter prior
s0_beta <- 2.75 # [1,inf); larger values = tighter prior
s100_alpha <- 10 # [1,inf]; larger values = tighter prior
s100_beta <- 1 # (0,1]; smaller values = tighter prior

bind_rows(prior(a = s0_alpha, b = s0_beta, label = "s0"),
          prior(a = s100_alpha, b = s100_beta, label = "s100")) %>% 
  ggplot(aes(x, y, col = factor(label), linetype = factor(label))) +
  geom_line() +
  # scale_colour_grey() +
  ylim(0,10) +
  labs(x = expression(theta), y = expression(paste('p(',theta,')', sep = '')),
       col = expression(theta), lty = expression(theta)) +
  theme(legend.position = c(0.2, 0.8))

# Model ----

setwd("~/great_escape/code")

# Ring sizes
ring_sizes <- 2.54 * c(3.5, 3.75, 4)

# Number of pots per set that were sampled for each treatment
ntrt_df <- counts %>% 
  group_by(Treatment, set) %>% 
  dplyr::summarise(ntrt = length(which(disposition == "tagged"))) %>% 
  data.table::dcast(set ~ Treatment)

# TMB data structure
data <- list(nset = length(unique(df$effort_no)), # number of sets
             nlen = length(unique(df$length_bin)), # number of length bins
             ntrt = length(unique(df$Treatment)), # number of treatments
             len = len, # vector of lengths for which there are data
             fit_len = seq(10, 100, 1), #fit_len, # vector of lengths for fitted values
             npot_ctl = ntrt_df$Control, # vector of number of control pots in each set
             npot_trt = as.matrix(select(ntrt_df, matches(paste(unique(df$Treatment), collapse = "|")))), # matrix of number of experimental pots in each set [j,k]
             ctl_dat = matrix(df$ctl_n[df$Treatment == "8.9 cm"], ncol = length(unique(df$effort_no))), #1), # # control pots: matrix of number of fish caught by length bin (row) by set (col)
             trt_dat = array(df$exp_n, # experimental pots: matrix of number of fish caught by length bin (row) by set (col)
                             dim = c(length(unique(df$length_bin)),
                                     length(unique(df$effort_no)),
                                     length(unique(df$Treatment)))),
             prior_type = 0, # 0 = normal, 1 = beta
             s0_index = s0_index, # Index of length where theoretical curves were 0 and 1 for each treatment
             s100_index = s100_index,
             sigma_s0 = 0.35,  # priors to constrain selectivity at 0 and 1 0.35
             sigma_s100 = 0.05, # 0.05
             s0_alpha = s0_alpha,
             s0_beta = s0_beta,
             s100_alpha = s100_alpha,
             s100_beta = s100_beta
) 

# Model parameters
parameters <- list(dummy = 0,
                   s50 = s50_vec,
                   slp = slp_vec,
                   log_delta = log(1),
                   nu = rep(0, length(unique(df$effort_no))))

# Bounds aka uniform priors in same order as parameter list
lb <- c(rep(0, 3), rep(0, 3), log(0.1), rep(-5, length(unique(df$effort_no))))
ub <- c(rep(100, 3), rep(0.9, 3), log(2), rep(5, length(unique(df$effort_no))))

# Compile
compile("escape.cpp")
dyn.load(dynlib("escape"))

# Model 1a ----

# Assume fixed delta with normal prior on s0 and s100
map <- list(dummy = factor(NA),
            log_delta = factor(NA))
data$prior_type <- 0 # normal
obj <- MakeADFun(data, parameters, map = map,
                   DLL = "escape", silent = TRUE,
                   hessian = TRUE, random = "nu")
opt <- nlminb(obj$par, obj$fn, obj$gr)

(Mod1a_AIC <- TMBAIC(opt))
best <- obj$env$last.par.best
(rep <- sdreport(obj))

plot_mle(delta = "fixed", selprior = "normal")

# Model 1b ----

# Estimate delta with normal prior on s0 and s0
map <- list(dummy = factor(NA))
data$prior_type <- 0 # normal
obj <- MakeADFun(data, parameters, map = map, 
                 DLL = "escape", silent = TRUE,
                 hessian = TRUE, random = "nu") 
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lb, upper = ub)

(Mod1b_AIC <- TMBAIC(opt))
best <- obj$env$last.par.best
(rep <- sdreport(obj))

plot_mle(delta = "estimated", selprior = "normal")

# Model 2a ----

# Assume fixed delta with beta prior on s0 and s100
map <- list(dummy = factor(NA),
            log_delta = factor(NA))
data$prior_type <- 1 # beta
obj <- MakeADFun(data, parameters, map = map,
                 DLL = "escape", silent = TRUE,
                 hessian = TRUE, random = "nu")
opt <- nlminb(obj$par, obj$fn, obj$gr)

(Mod2a_AIC <- TMBAIC(opt))
best <- obj$env$last.par.best
(rep <- sdreport(obj))
Mod2a_nll <- obj$fn()[1]
plot_mle(delta = "fixed", selprior = "beta")

# Model 2b ----

# Estimate delta with beta prior on s0 and s0
map <- list(dummy = factor(NA))
data$prior_type <- 1 # beta
obj <- MakeADFun(data, parameters, map = map, 
                 DLL = "escape", silent = TRUE,
                 hessian = TRUE, random = "nu") 
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lb, upper = ub)

(Mod2b_AIC <- TMBAIC(opt))
best <- obj$env$last.par.best
(rep <- sdreport(obj))
Mod2b_nll <- obj$fn()[1]

plot_mle(delta = "estimated", selprior = "beta")

# AIC model selection ----
AICtable <- data.frame(Model = c("Mod1a_fixdelta_normal",
                                 "Mod1b_estdetla_normal",
                                 "Mod2a_fixdelta_beta",
                                 "Mod2b_estdetla_beta"),
                       AIC = c(Mod1a_AIC, Mod1b_AIC, Mod2a_AIC, Mod2b_AIC))
AICtable
# Mod 2a and 2b are the same by AIC, use one that has lowest nll
Mod2a_nll; Mod2b_nll

# MCMC ----

library(tmbstan)

# Run in parallel with a init function
cores <- parallel::detectCores()-1
options(mc.cores = cores)
init.fn <- function(){
  list(s50 = sort(runif(3, 55, 70)), 
       slp = sort(runif(3, 0.1, 0.2), decreasing = TRUE),
       log_delta = log(rnorm(1, 1, 0.1)),
       nu = rnorm(length(unique(df$effort_no))))}

fit <- tmbstan(obj, seed = 1, chains = 3, 
               iter = 10000,
               warmup = 2000, thin = 10,
               open_progress = FALSE, 
               init = init.fn, lower = lb, upper = ub)

pdf(file = "../figures/pairs.pdf", width = 7.08, height = 7.08/1.618)#, width = 180, height = 180)# , dpi = 600, units = "mm")
pairs(fit, pars = names(obj$par)) # Pairs plot of the fixed effects
dev.off()

# Explore the fit use shinystan
# launch_shinystan(fit)

## Can also get ESS and Rhat from rstan::monitor
mon <- monitor(fit)
write_csv(mon, "../output/convergence.csv")
max(mon$Rhat)
min(mon$Bulk_ESS)
min(mon$Tail_ESS)

# Other methods provided by 'rstan'
class(fit)
methods(class="stanfit")

# Trace plot
dev.new()
trace <- traceplot(fit, pars = names(obj$par), inc_warmup = TRUE)
trace + scale_color_grey() + theme(legend.position = c(0.7,0.15), legend.direction = "horizontal")
ggsave(filename = "../figures/trace.pdf", width = 180, height = 180/1.618, units = "mm")

# Extract marginal posteriors easily
post <- as.matrix(fit)

# hist(post[,'nu[1]'])                     # random effect
# hist(post[,'s50[1]'])                    # fixed effect
# dim(post)
write_csv(as.data.frame(post), "../output/posterior_samples.csv")

# Summary of parameter estimates 
pars_sum <- summary(fit)$summary
write_csv(as.data.frame(pars_sum), "../output/param_summary.csv")

# Posterior for derived quantities slx and phi. The last column in post is the
# log-posterior density (lp__) and needs to be dropped via -ncol(post)
slx <- list()
phi <- list()

for(i in 1:nrow(post)){
  r <- obj$report(post[i,-ncol(post)])
  slx[[i]] <- cbind(r$full_slx, rep(i, nrow(r$full_slx)), fit_len)
  phi[[i]] <- cbind(r$fit_phi, rep(i, nrow(r$fit_phi)), len)
}

slx <- as.data.frame(do.call(rbind, slx))
names(slx) <- c(paste(unique(df$Treatment)), "iter", "length_bin")
slx <- slx %>% 
  pivot_longer(-c(iter, length_bin), names_to = "Treatment", values_to = "slx") %>% 
  group_by(Treatment, length_bin) %>% 
  summarize(mean = mean(slx),
            median = median(slx),
            q025 = quantile(slx, 0.025),
            q975 = quantile(slx, 0.975))

phi <- as.data.frame(do.call(rbind, phi))
names(phi) <- c(paste(unique(df$Treatment)), "iter", "length_bin")
phi <- phi %>% 
  pivot_longer(-c(iter, length_bin), names_to = "Treatment", values_to = "phi") %>% 
  group_by(Treatment, length_bin) %>% 
  summarize(mean = mean(phi),
            median = median(phi),
            q025 = quantile(phi, 0.025),
            q975 = quantile(phi, 0.975))

write_csv(slx, "../output/slx_ci.csv")
write_csv(phi, "../output/phi_ci.csv")

com %>% 
  select(Treatment, length_bin, p = combined_p) %>% 
  left_join(phi, by = c("Treatment", "length_bin")) %>% 
  mutate(resid = p - mean,
         Treatment = factor(Treatment, levels = c("8.9 cm", "9.5 cm", "10.2 cm"), 
                            ordered = TRUE)) -> phi

# Figures ----

slx %>% 
  ungroup() %>%
  mutate(Method = "SELECT") %>% 
  bind_rows(sprior %>% 
              rename(mean = p, length_bin = length) %>% 
              mutate(Method = "Theoretical (May/Jun)") %>% 
              filter(Treatment != "Control" &
                       length_bin %in% seq(30, 100, 0.4))) %>% 
  mutate(Treatment = factor(Treatment, levels = c("8.9 cm", "9.5 cm", "10.2 cm"), 
                            ordered = TRUE),
         Method = factor(Method, levels = c("Theoretical (May/Jun)", "SELECT"),
                         ordered = TRUE)) -> slx

slx %>% filter(length_bin == 63)
slx %>% 
  filter(q025 > 0.99 & Method == "SELECT") %>% 
  group_by(Treatment) %>% 
  filter(length_bin == min(length_bin))

# Selectivity
p_slx <- ggplot() +
  geom_vline(xintercept = 63, col = "grey90", lty = 5, size = 0.2) + #
  geom_line(data = slx, aes(x = length_bin, y = mean, linetype = Treatment, 
                            colour = Method, size = Method,
                            group = interaction(Method, Treatment))) +
  geom_ribbon(data = slx %>% 
                filter(Method == "SELECT"),
              aes(x = length_bin, group = Treatment, ymin = q025, ymax = q975), 
              col = NA, alpha = 0.1, show.legend = FALSE) +
  scale_colour_manual(values = c("grey60", "grey10")) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  scale_size_manual(values = c(0.4, 0.6)) +
  # annotate("curve", x = 45, y = 0.85, xend = 62.5, yend = 1, size = 0.2,
  #          colour = "grey70", curvature = -0.3, arrow = arrow(length = unit(1, "mm"))) +
  # annotate("text", x = 45, y = 0.8, colour = "grey60", size = 3,
  #          label = as.character(expression(paste(italic(L)[50]== "63 cm"))), parse = TRUE) +
  xlim(30, 90) +
  labs(x = "Fork length (cm)", y = "Proportion retained") + 
  theme(#legend.position = c(0.8, 0.2),
        # legend.key.width = unit(1.6,"line"),
        legend.position = c(0.79, 0.3),
        # legend.text=element_text(size = 7),
        legend.spacing.y = unit(0, "cm"))
p_slx
ggsave(filename = "../figures/slx_ci.pdf", plot = p_slx, device = "pdf",
       dpi = 600, units = "mm", width = 80, height = 80/1.618)

# Phi residuals
p_resid <- ggplot(phi, aes(x = length_bin, y = resid)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_segment(aes(x = length_bin, xend = length_bin, y = 0, yend = resid),
               size = 0.2, col = "grey") +
  geom_point() +
  facet_wrap(~Treatment, ncol = 1) +
  labs(x = "Fork length (cm)", y = "Residuals") +
  theme(strip.text.x = element_text(size=0))

# Phi 
phi <- phi %>% mutate(Treatment2 = Treatment)
label_df <- distinct(phi, Treatment)
p_phi <- ggplot(phi, aes(x = length_bin, group = Treatment)) + 
    geom_line(data = select(phi, -Treatment),
            aes(y = mean, group = Treatment2), col = "grey") +
  geom_ribbon(aes(ymin = q025, ymax = q975),
              alpha = 0.1, col = NA) +
  geom_point(aes(y = p)) + 
  geom_line(aes(y = mean)) + 
  labs(x = "Fork length (cm)", y = "Proportion in treatment pots") +
  geom_text(data = label_df, size = 2.5, fontface = "plain", aes(x = 55, y = 0.6, label = Treatment)) +
  facet_wrap(~ Treatment, ncol = 1) +
  ylim(c(0,0.7)) +
  theme(strip.text.x = element_text(size=0))

p <- plot_grid(p_phi, p_resid, ncol = 2)
p
ggsave(filename = "../figures/phi_resids.pdf", plot = p, dpi = 600, 
       device = "pdf", units = "mm", height = 180, width = 180)

# Generalized curves -----

# Regress s50 and k estimates on escape ring diameter to develop generalized
# selectivity equation (Eq 6 in Arana et al 2011) *Note that our
# parameterization of the logistic model is slightly different.

s50_est <- pars_sum[1:3,1]
slp_est <- pars_sum[4:6,1]

s50_lm <- lm(s50_est ~ ring_sizes)
coef(s50_lm)
summary(s50_lm)

slp_lm <- lm(slp_est ~ ring_sizes)
coef(slp_lm)
summary(slp_lm)

ring_vec <- seq(8.5, 10.5, 0.4)

p_gen <- matrix(nrow = length(fit_len), ncol = length(ring_vec))

for(i in 1:length(fit_len)) {
  for(j in 1:length(ring_vec)) {
  p_gen[i,j] <- 1 / (1 + exp(-1 * (coef(slp_lm)[1] + coef(slp_lm)[2] * ring_vec[j]) *
                               (fit_len[i] - (coef(s50_lm)[1] + coef(s50_lm)[2] * ring_vec[j]))))
  }
}

p_gen <- as.data.frame(p_gen)
names(p_gen) <- paste(ring_vec, "cm")
p_gen <- p_gen %>% 
  mutate(length_bin = fit_len) %>% 
  pivot_longer(-length_bin, names_to = "Escape ring", values_to = "p")
  
p <- ggplot(p_gen, aes(x = length_bin, y = p, col = `Escape ring`, group = `Escape ring`)) +
  geom_line() +
  scale_color_grey() +
  labs(x = "Fork length (cm)", y = "Proportion retained") +
  xlim(c(35, 80)) +
  theme(legend.position = c(0.79, 0.3),
        legend.text=element_text(size = 6),
        legend.spacing.y = unit(0, "cm"))
p
ggsave(filename = "../figures/gen_slx.pdf", plot = p, device = "pdf",
       dpi = 600, units = "mm", width = 80, height = 80/1.618)

# Evaluation of beta priors ----

prior_0 
prior_100

# Comparison of priors and posterior samples on 0 and 100% retention probs

post_priors <- list()

for(i in 1:nrow(post)){
  r <- obj$report(post[i,-ncol(post)])
  # Extract selectivity probability (proportion retained) at each of the mean
  # prior lengths
  post_priors[[i]] <- tibble(cbind(r$fit_slx, rep(i, nrow(r$fit_slx)), tmb_index)) %>% 
    filter(index %in% s0_index) %>% 
    mutate(prior = "0") %>% 
    bind_rows(tibble(cbind(r$fit_slx, rep(i, nrow(r$fit_slx)), tmb_index)) %>% 
                filter(index %in% s100_index) %>% 
                mutate(prior = "1")
    )
  }

post_priors <- as.data.frame(do.call(rbind, post_priors))
names(post_priors) <- c(paste(unique(df$Treatment)), "iter", "length_bin", "index", "prior")
post_priors <- post_priors %>% 
  pivot_longer(-c(iter, length_bin, index, prior), names_to = "Treatment", values_to = "p")# %>% 

prior_lkup <- prior_0 %>% 
  mutate(prior = "0") %>% 
  bind_rows(prior_100 %>% 
              mutate(prior = "1")) %>% 
  mutate(id = paste(Treatment, length_bin, prior, sep = "_")) 

post_priors <- post_priors %>% 
  mutate(id = paste(Treatment, length_bin, prior, sep = "_")) %>% 
  filter(id %in% prior_lkup$id)

ggplot() +
  geom_histogram(data = post_priors, aes(x = p, col = prior, fill = prior),
                 stat = "density") +
  facet_wrap(~Treatment) +
  geom_line(data = bind_rows(prior(a = s0_alpha, b = s0_beta, label = "0"),
                             prior(a = s100_alpha, b = s100_beta, label = "1")) %>% 
              rename(prior = label),
            aes(x, y, col = prior, linetype = prior))
