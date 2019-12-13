# Escape ring selectivity model
# Jane Sullivan
# jane.sullivan1@alaska.gov
# Last updated 2019-11-26

source("code/helper.r")
library(TMB)
library(cowplot)
# library(nimble) # inverse gamma prior
library(broom) # tidy() clean up parameter estimates
library(tmbstan) # run mcmc
library(shinystan) # awesome diagnostic tool

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
# the control pot is 100% selected
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
  ylim(c(0,1)) +
  geom_smooth()

ggplot(com, aes(x = length_bin, y = combined_p,
                   group =  Treatment, col = Treatment)) +
  geom_point() +
  ylim(c(0,1)) +
  geom_smooth()

# Priors ----
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
# s100_index <- c(17,22,26)

# Get starting values from theoretical selectivity curves
s50_vec <- vector(length = 3)
slp_vec <- vector(length = 3)

for(i in 1:length(unique(df$Treatment))) {
  tmp <- sprior %>% filter(Treatment == unique(df$Treatment)[i])
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

s0_alpha <- 0.7
s0_beta <- 20
s100_alpha <- 40.0
s100_beta <- 0.1

bind_rows(prior(a = s0_alpha, b = s0_beta, label = "s0"),
          prior(a = s100_alpha, b = s100_beta, label = "s100")) %>% 
  ggplot(aes(x, y, col = factor(label), linetype = factor(label))) +
  geom_line(size = 1) +
  scale_colour_grey() +
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
  dcast(set ~ Treatment)

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

# Bounds aka uniform priors
lb <- c(rep(0, 3), rep(0, 3), log(0.1), rep(-5, length(unique(df$effort_no))))
ub <- c(rep(100, 3), rep(0.9, 3), log(2), rep(5, length(unique(df$effort_no))))

# Compile
compile("escape.cpp")
dyn.load(dynlib("escape"))

# Model 1: Assume fixed delta
# map <- list(dummy = factor(NA),
#             log_delta = factor(NA))
# 
# obj <- MakeADFun(data, parameters, map = map, 
#                    DLL = "escape", silent = TRUE,
#                    hessian = TRUE, random = "nu") 
# 
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# 
# Mod1_AIC <- TMBAIC(opt)

# Model 2: Estimate delta
map <- list(dummy = factor(NA))

obj <- MakeADFun(data, parameters, map = map, 
                 DLL = "escape", silent = TRUE,
                 hessian = TRUE, random = "nu") 

opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lb, upper = ub)

# Mod2_AIC <- TMBAIC(opt)
# 
# Mod1_AIC - Mod2_AIC

best <- obj$env$last.par.best
print(best)
rep <- sdreport(obj)
print(rep)

obj$report()$prior_s0
obj$report()$prior_s100
obj$report()$set_effect

phi <- as.data.frame(obj$report()$fit_phi)
names(phi) <- paste0(unique(df$Treatment))
phi <- phi %>% mutate(length_bin = len)
phi <- melt(data = phi, id.vars = "length_bin", variable.name = "Treatment", value.name = "phi") 

com %>% 
  left_join(phi, by = c("Treatment", "length_bin")) %>% 
  mutate(resid = combined_p - phi,
         Treatment = factor(Treatment, levels = c("8.9 cm", "9.5 cm", "10.2 cm"), 
                            ordered = TRUE)) -> tmp

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

slx <- as.data.frame(obj$report()$full_slx)
names(slx) <- paste0(unique(df$Treatment))
slx <- slx %>% mutate(length_bin = fit_len)# sort(unique(df$length_bin)))
slx <- melt(data = slx, id.vars = "length_bin", variable.name = "Treatment", value.name = "slx")

ggplot() +
  ylim(c(0,1)) +
  geom_line(data = slx, aes(x = length_bin, y = slx, group = Treatment, col = Treatment))

# MCMC ----

# fit <- tmbstan(obj, chains = 1)

# Run in parallel with a init function
cores <- parallel::detectCores()-1
options(mc.cores = cores)
init.fn <- function(){
  list(s50 = sort(runif(3, 55, 70)), slp = sort(runif(3, 0.1, 0.2), decreasing = TRUE),
       log_delta = log(rnorm(1, 1, 0.1)),
       nu = rnorm(length(unique(df$effort_no))))}

fit <- tmbstan(obj, chains = cores, open_progress = FALSE, init = init.fn, lower = lb, upper = ub)

pdf(file = "../figures/pairs.pdf", width = 7.08, height = 7.08)#, width = 180, height = 180)# , dpi = 600, units = "mm")
pairs(fit, pars = names(obj$par)) # Pairs plot of the fixed effects
dev.off()

# Explore the fit use shinystan
launch_shinystan(fit)

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
ggsave(filename = "../figures/trace.pdf", width = 180, height = 180, units = "mm")

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
  melt(id.vars = c("iter", "length_bin"), variable.name = "Treatment", value.name = "slx") %>% 
  group_by(Treatment, length_bin) %>% 
  summarize(mean = mean(slx),
            median = median(slx),
            q025 = quantile(slx, 0.025),
            q975 = quantile(slx, 0.975))

phi <- as.data.frame(do.call(rbind, phi))
names(phi) <- c(paste(unique(df$Treatment)), "iter", "length_bin")
phi <- phi %>% 
  melt(id.vars = c("iter", "length_bin"), variable.name = "Treatment", value.name = "phi") %>% 
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
              mutate(Method = "Theoretical (May)") %>% 
              filter(Treatment != "Control" &
                       length_bin %in% seq(30, 100, 0.4))) %>% 
  mutate(Treatment = factor(Treatment, levels = c("8.9 cm", "9.5 cm", "10.2 cm"), 
                            ordered = TRUE),
         Method = factor(Method, levels = c("Theoretical (May)", "SELECT"),
                         ordered = TRUE)) -> slx

# Selectivity
p_slx <- ggplot() +
  geom_vline(xintercept = 61, col = "grey90", lty = 5, size = 0.2) + #
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
  annotate("curve", x = 45, y = 0.85, xend = 60.5, yend = 0.99, size = 0.2,
           colour = "grey70", curvature = -0.3, arrow = arrow(length = unit(1, "mm"))) +
  annotate("text", x = 45, y = 0.8, colour = "grey60", size = 3,
           label = as.character(expression(paste(italic(L)[50]== "61 cm"))), parse = TRUE) +
  xlim(30, 90) +
  labs(x = "Fork length (cm)", y = "Proportion retained") + 
  theme(#legend.position = c(0.8, 0.2),
        # legend.key.width = unit(1.6,"line"),
        legend.position = c(0.79, 0.3),
        legend.text=element_text(size = 7),
        legend.spacing.y = unit(0, "cm"))
p_slx
ggsave(filename = "../figures/slx_ci.pdf", plot = p_slx, device = "pdf",
       dpi = 600, units = "mm", width = 80, height = 80)

# Phi residuals
p_resid <- ggplot(phi, aes(x = length_bin, y = resid)) +
  geom_hline(yintercept = 0, col = "grey", size = 1) +
  geom_segment(aes(x = length_bin, xend = length_bin, y = 0, yend = resid),
               size = 0.2, col = "grey") +
  geom_point() +
  facet_wrap(~Treatment, ncol = 1) +
  labs(x = "Fork length (cm)", y = "Residuals") +
  theme(strip.text.x = element_text(size=0))

slx %>% filter(Method == "SELECT" & length_bin == 61.0)

# Phi 
phi <- phi %>% mutate(Treatment2 = Treatment)
p_phi <- ggplot(phi, aes(x = length_bin, group = Treatment)) + 
    geom_line(data = select(phi, -Treatment),
            aes(y = mean, group = Treatment2), col = "grey") +
  geom_ribbon(aes(ymin = q025, ymax = q975),
              alpha = 0.1, col = NA) +
  geom_point(aes(y = p)) + 
  geom_line(aes(y = mean)) + 
  labs(x = "Fork length (cm)", y = "Proportion in treatment pots") +
  geom_text(aes(x = 55, y = 0.6, label = Treatment)) +
  facet_wrap(~ Treatment, ncol = 1) +
  theme(strip.text.x = element_text(size=0))

p <- plot_grid(p_phi, p_resid, ncol = 2)
p
ggsave(filename = "../figures/phi_resids.pdf", plot = p, dpi = 600, 
       device = "pdf", units = "mm", height = 150, width = 180)

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
  melt(id.vars = "length_bin", variable.name = "Escape ring", value.name = "p")

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
       dpi = 600, units = "mm", width = 80, height = 80)
