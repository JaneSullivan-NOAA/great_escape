# Simple logistic regression
# Jane Sullivan
# jane.sullivan1@alaska.gov
# Last updated 2019-10-07

source("code/helper.r")
library(TMB)
library(cowplot)

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

# Length bins
# Same length bin structure used by Feds)
bio <- bio %>%
  filter(!c(length < 40)) %>%
  mutate(length2 = ifelse(length < 41, 41,
                          ifelse(length > 99, 99, length)),
         length_bin = cut(length2, breaks = seq(39.9, 99.9, 2),
                          labels = paste(seq(41, 99, 2)))) %>%
  select(-length2)

# bio <- bio %>%
#   filter(!c(length < 40)) %>%
#   mutate(length2 = ifelse(length < 41, 41,
#                           ifelse(length > 79, 79, length)),
#          Length_bin = cut(length2, breaks = seq(39.9, 79.9, 2),
#                           labels = paste(seq(41, 79, 2)))) %>%
#   select(-length2)

# Reorganize data and get number of fish caught in experimental pots / control
# pots for each length bin
sum_df <- bio %>% 
  group_by(Treatment, effort_no, length_bin, .drop=FALSE) %>%  #
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
         # FLAG for exp_n / ctl_n:  clean up the data here. If ctl pot = 0, make
         # p2 = exp_n / ctl_n,
         # # p = 0. If exp > ctl, make p = 1
         # p = ifelse(p2 > 1, 1, ifelse(is.na(p2), 0, p2))

ggplot(sum_df, aes(x = length_bin, y = p,
                   group =  Treatment, col = Treatment)) +
  geom_point() +
  geom_line()

ggplot(sum_df, aes(x = length_bin, y = p,
                   group =  factor(effort_no), col = factor(effort_no))) +
  geom_point() +
  geom_line() +
  facet_wrap(~Treatment)

TREATMENT <- "4.00 in"
# TREATMENT <- "3.50 in"
trt <- filter(sum_df, Treatment == TREATMENT)

# trt <- mutate(trt, Effort_no = factor(effort_no))
# trt <- mutate(trt, effort_no = 1)

# mod <- glm(p ~ length_bin, family = "binomial", weights = tot_n, data = trt)
# mod <- glmer(p ~ length_bin + (1 | effort_no), family = "binomial", weights = tot_n, data = trt)
# tidy(mod)
# summary(mod)
# pred_df <- data.frame(length_bin = seq(30, 100, 1)) %>% 
#   mutate(effort_no = 6)
# pred_df$pred <- predict(mod, pred_df, type = "response")

# Data, parameter, and map for TMB ----

data <- list(model = 1, # model switch
             nset = length(unique(trt$effort_no)), # number of sets
             nlen = length(unique(trt$length_bin)), # number of length bins
             len = sort(unique(trt$length_bin)), # vector of lengths
             npot_ctl = rep(4, length(unique(trt$effort_no))), # vector of number of control pots in each set
             npot_exp = rep(4, length(unique(trt$effort_no))), # vector of number of experimental pots in each set
             ctl_dat = matrix(trt$ctl_n, ncol = length(unique(trt$effort_no))), # control pots: matrix of number of fish caught by length bin (row) by set (col)
             exp_dat = matrix(trt$exp_n, ncol = length(unique(trt$effort_no)))) # experimental pots: matrix of number of fish caught by length bin (row) by set (col)

# Starting values from Table N.4 Model 1.3 and 1.2 in Haist et al. 2004
parameters <- list(dummy = 0,
                   log_s50 = log(60),
                   # log_slxr = log(10),
                   # log_lslx = log(64-50),
                   # log_uslx = log(69-64),
                   log_s90 = log(69),
                   log_s10 = log(50),
                   log_delta = log(1),
                   nu = rep(0, length(unique(trt$effort_no))))

# Troubleshooting map
map <- list(log_s50 = factor(NA),
            log_slxr = factor(NA),
            # log_lslx = factor(NA),
            # log_uslx = factor(NA),
            log_delta = factor(NA),
            nu = rep(factor(NA), length(unique(trt$effort_no))))

# Compile TMB code and fit model ----
setwd("~/great_escape/code")

compile("escape.cpp")
dyn.load(dynlib("escape"))

map <- list(dummy = factor(NA),
            # log_s50 = factor(NA),
            log_lslx = factor(NA)
            # log_uslx = factor(NA),
            # log_delta = factor(NA),
            # nu = rep(factor(NA), length(unique(trt$effort_no)))
            )

lowbnd= c(log(45), # log_s50
          # log(2),  # log_slxr
          # log(1), # log_lslx
          log(1), # log_uslx
          log(0.5), # log_delta
          rep(-5, data$nset)
          ) 

uppbnd= c(log(75), # log_s50
          # log(25),  # log_slxr
          # log(20), # log_lslx
          log(20), # log_uslx
          log(1.5),  # log_delta
          rep(5, data$nset)
          )  

data$model <- 1
model <- MakeADFun(data, parameters, map = map, 
                   DLL = "escape", silent = FALSE,
                   hessian = TRUE, random = "nu") #
# checking for minimization
xx <- model$fn(model$env$last.par)
print(model$report())

fit <- nlminb(model$par, model$fn, model$gr, 
              # control=list(eval.max=1000000,iter.max=100000),
              control=list(rel.tol=1e-12,
                           eval.max=100000,iter.max=10000),
              lower=lowbnd,upper=uppbnd)
best <- model$env$last.par.best
print(best)
exp(best[3])
rep <- sdreport(model)
print(rep)

s50 <- model$report()$s50
s10 <- model$report()$s10
s90 <- model$report()$s90
# r <- exp(best[2])# model$report()$slxr
phi <- model$report()$phi

slx <- vector(length = data$nlen)

for(i in 1:data$nlen) {
  len <- data$len[i]
  if(len <= s50) {
    slx[i] <- 1 / (1 + exp(-2 * log(3) * ((len - s50) / (s50 - s10))))
  } else {
    slx[i] <- 1 / (1 + exp(-2 * log(3) * ((len - s50) / (s90 - s50))))
  }
}

plot(data$len, slx)


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
