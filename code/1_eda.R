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

# Fishery girth data (collected in order to account for grith seasonal adjustments)
fsh_grth <- read_csv(paste0("data/fsh_bio_", YEAR, ".csv")) 

# Effort data 
effort <- read_csv(paste0("data/pot_effort_", YEAR, ".csv"))

# Total sablefish counts by depth and disposition (some pots in a set were
# dumped due to processing time)
counts <- read_csv(paste0("data/total_counts_", YEAR, ".csv")) %>% 
  mutate(Treatment = derivedFactor("Control" = treatment == "Blue",
                                   "3.50 in" = treatment == "Purple",
                                   "3.75 in" = treatment == "Green",
                                   "4.00 in" = treatment == "Yellow",
                                   .default = NA,
                                   .ordered = TRUE))

# Email sent to A. Baldwin 2019-10-01 about the NAs. See GitHub issue #1
bio <- bio %>% filter(!is.na(Treatment))

# Length frequency ----

p <- bio %>% 
  filter(!is.na(Treatment) & between(length, 40, 90)) %>%
  droplevels() %>%
  ggplot(aes(x = length, colour = Treatment, size = Treatment, linetype = Treatment)) + 
  geom_freqpoly() +
  scale_colour_manual(values = c("grey90", "grey70", "grey40", "black")) +
  scale_size_manual(values = c(1.5, 0.7, 0.7, 0.7)) +
  scale_linetype_manual(values = c(1, 1, 2, 3)) +
  xlim(40, 90) +
  labs(x = "\nLength (cm)", y = "Count\n") + 
  theme(legend.position = c(0.8, 0.7))

# Caption: Length frequency distribution by escape ring treatment.
ggsave(plot = p, filename = paste0("figures/size_freq_", YEAR, ".png"), dpi=300, height=3, width=6, units="in")

# Girth outliers  ----

# Identify outliers, plot them, remove them
grth <- bio %>% filter(!is.na(girth))
# identify(length, girth, plot=TRUE)

grth <- grth %>% 
  mutate(Outlier = ifelse(row_number() %in% c(190, 532, 533), "Outlier", "Normal"))

p <- ggplot(grth, aes(x = length, y = girth, col = Outlier, shape = Outlier)) +
  geom_point() +
  scale_colour_grey() +
  scale_shape_manual(values = c(20, 8)) +
  labs(x = "\nLength (cm)", y = "Girth (mm)\n") +
  theme(legend.position = c(0.8, 0.2))

ggsave(plot = p, filename = paste0("figures/girth_outliers_", YEAR, ".png"), 
       dpi=300, height=3, width=6, units="in")

grth <- grth %>% filter(Outlier != "Outlier") %>% 
  select(-Outlier)

# Plot cleaned up girth measurements by treatment
p <- ggplot(grth, aes(x = length, y = girth)) +
  geom_point(shape = 20) +
  facet_wrap(~ Treatment) +
  labs(x = "\nLength (cm)", y = "Girth (mm)\n")

ggsave(plot = p, filename = paste0("figures/girth_bytreatment_", YEAR, ".png"), 
       dpi=300, height=5, width=6, units="in")

# Girth by treatment ----

# Girth isn't normally distributed. Log-transformation improves it but
# assumption still violated. Consequently I used I used glm instead of lm
ggdensity(log(grth$girth))
ggqqplot(log(grth$girth))
shapiro.test(log(grth$girth))

# Fit glm first with all treatments combined, then with treatments separate.
# Log-transform girth and length. Use AIC to determine necessary complexity
fit_simple <- glm(log(girth) ~ log(length), 
                  family = gaussian(link = "identity"), data = grth)
fit_complex <- glm(log(girth) ~ log(length) * factor(Treatment, ordered = FALSE), 
                   family = gaussian(link = "identity"), data = grth)
AIC(fit_simple, fit_complex) # more complex model was favored by a delta AIC of 3.5

# Model summaries and goodness of fit R2
summary(fit_simple)
summary(fit_complex)
R2_simple <- 1 - (fit_simple$deviance / fit_simple$null.deviance)
R2_complex <- 1 - (fit_complex$deviance / fit_complex$null.deviance)

# Get fitted values for girth and prediction intervals (need to exponentiate for figure)
pred <- bio %>% select(Treatment, length)
pred <- ciTools::add_pi(pred, fit_complex, alpha = 0.05, names = c("pi_lwr", "pi_upp"))

# Apply bias correction
pred <- pred %>% 
  mutate(fitted = exp(pred) * exp(0.5 * sigma(fit_complex)^2),
         lower = exp(pi_lwr) * exp(0.5 * sigma(fit_complex)^2),
         upper = exp(pi_upp) * exp(0.5 * sigma(fit_complex)^2))

p <- ggplot() +
  geom_point(data = grth, aes(x = length, y = girth), shape = 20) +
  geom_ribbon(data = pred, aes(x = length, ymin = lower, ymax = upper), 
              alpha = 0.6, fill = "grey70") +
  geom_line(data = pred, aes(x = length, y = fitted, group = Treatment)) +
  facet_wrap(~ Treatment) +
  labs(x = "\nLength (cm)", y = "Girth (mm)\n") 
p
ggsave(plot = p, filename = paste0("figures/fitted_girth_bytreatment_", YEAR, ".png"), 
       dpi=300, height=5, width=6, units="in")

# Girth adjustments  ----

# Combine girths from survey and fishery
comb_grth <- grth %>% 
  filter(Treatment == "Control") %>%
  select(length, girth) %>% 
  mutate(Source = "Survey (May)") %>% 
  bind_rows(fsh_grth %>% 
              select(length, girth) %>% 
              mutate(Source = "Fishery (Sep and Oct)"))

fit_simple <- glm(log(girth) ~ log(length), 
                  family = gaussian(link = "identity"), data = comb_grth)
fit_int <- glm(log(girth) ~ log(length) + Source, 
               family = gaussian(link = "identity"), data = comb_grth)
fit_intslp <- glm(log(girth) ~ log(length) * Source, 
                  family = gaussian(link = "identity"), data = comb_grth)
AIC(fit_simple, fit_int, fit_intslp) # best model = different intercepts

summary(fit_int)
R2_int <- 1 - (fit_int$deviance / fit_int$null.deviance)

# Get fitted values for girth and prediction intervals (need to exponentiate for figure)
pred <- comb_grth %>% select(Source, length)
pred <- ciTools::add_pi(pred, fit_int, alpha = 0.05, names = c("pi_lwr", "pi_upp"))

# Apply bias correction
pred <- pred %>% 
  mutate(fitted = exp(pred) * exp(0.5 * sigma(fit_int)^2),
         lower = exp(pi_lwr) * exp(0.5 * sigma(fit_int)^2),
         upper = exp(pi_upp) * exp(0.5 * sigma(fit_int)^2))

p <- ggplot() +
  geom_ribbon(data = pred, aes(x = length, ymin = lower, ymax = upper, fill = Source), 
              alpha = 0.3) +
  geom_point(data = comb_grth, aes(x = length, y = girth, colour = Source, shape = Source), size = 0.8) +
  geom_line(data = pred, aes(x = length, y = fitted, group = Source, colour = Source, linetype = Source), size = 1) +
  scale_colour_manual(values = c("grey10", "grey60")) +
  scale_fill_manual(values = c("grey80", "grey70")) +
  labs(x = "\nLength (cm)", y = "Girth (mm)\n") +
  theme(legend.position = c(0.8, 0.2))

# Caption: A comparison of fitted values and prediction intervals for the
# regression of girth on length for data collected during the survey in May
# (grey triangles) and fishery in September and October (black circles).
ggsave(plot = p, filename = paste0("figures/girth_bysource_", YEAR, ".png"), 
       dpi=300, height=3, width=6, units="in")

# Theoretical selectivity curves ----

# Approach similar to Treble et al. 1998 used for Southern rock lobster in
# Australia (Fisheries Research 34 1998 289–305)

# Potential development - add normal distribution curves to girth-length
# regression to illustrate method used to obtain theoretical selectivity curves:
# https://stackoverflow.com/questions/31794876/ggplot2-how-to-curve-small-gaussian-densities-on-a-regression-line

# For theoretical selectivity curves, use model output from Control data only.
grth <- grth %>% filter(Treatment == "Control")
fit <- glm(log(girth) ~ log(length), family = gaussian(link = "identity"), data = grth)
summary(fit)
girth_se <- sigma(fit) # se of girth

# fit <- glm(girth ~ length, family = gaussian(link = "identity"), data = grth)
# summary(fit)
# girth_se <- sigma(fit)

pred_df <- data.frame(length = seq(30, 100, 0.01))
pred_df$pred <- predict(fit, pred_df)

# A. Baldwin measured the mesh size diameter for me. Assume 73 mm. See issue 2
# on github for documentation.
ring <- data.frame(ring_in = c(73 / 25.4, 3.5, 3.75, 4)) %>% 
  mutate(ring_mm = ring_in * 25.4)

# Assume girths are lognormally distributed. Simulate girth distribution at 1 cm
# increments, divide by pi to get approximate fish diameter. Determine proportion
# retained for each treatment, assuming any fish diameter < the ring diameter
# would be able to escape.
sel <- matrix(nrow = length(pred_df$pred),
              ncol = length(ring$ring_mm))
nsim <- 5000

for(i in 1:length(pred_df$pred)) {
  
  sim <- exp(rnorm(n = nsim, mean = pred_df$pred[i], sd = girth_se)) / pi
  
  for(j in 1:length(ring$ring_mm)) {
    
    sel[i,j] <- length(which(sim > ring$ring_mm[j])) / nsim
  }
}

# Reformat and plot
sel <- as.data.frame(sel)
names(sel) <- levels(bio$Treatment)
sel <- sel %>% mutate(length = pred_df$length)
sel <- data.table::melt(data = sel, id.vars = "length", variable.name = "Treatment", value.name = "p")

# Save output, 90% selectivity used to inform prior
write_csv(sel, paste0("output/theoretical_selectivity_", YEAR, ".csv"))

sel <- sel %>% filter(length %in% seq(30, 100, 0.1))
p <- ggplot(sel, aes(x = length, y = p, col = Treatment, 
                     linetype = Treatment, group = Treatment, size = Treatment)) +
  geom_line() +
  scale_colour_manual(values = c("grey90", "grey70", "grey40", "black")) +
  scale_size_manual(values = c(1.5, 0.7, 0.7, 0.7)) +
  scale_linetype_manual(values = c(1, 1, 2, 3)) +
  labs(x = "\nLength (cm)", y = "Proportion retained\n") +
  theme(legend.position = c(0.8, 0.3))
p
# Caption: Theoretical selectivity curves for control and escape ring treatments.
ggsave(plot = p, filename = paste0("figures/theoretical_selectivity_", YEAR, ".png"),
       dpi=300, height=3, width=6, units="in")

# Theoretical w/ soak time ----

# Filter out control selectivity
ctl <- sel %>% filter(Treatment == "Control")

# Probability of a fish finding the escape ring as a function of soak time:
soak <- data.frame(hr = seq(1, 50, 1))
x50 <- 20
x95 <- 40
soak <- soak %>% 
  mutate(p = 1 / (1 + exp(-log(19) * (hr - x50) / (x95 - x50))))
p <- ggplot(soak, aes(x = hr, y = p)) +
  geom_line() +
  ylim(c(0, 1)) +
  labs(x = "\nSoak time (hr)",
       y = "Probability of finding escape ring\n")
p
# Caption: Theoretical probability of a fish finding an escape ring as a function of soak time.
ggsave(plot = p, filename = paste0("figures/escape_prob_soaktime.png"),
       dpi=300, height=3, width=6, units="in")

pred_df <- data.frame(length = seq(30, 100, 0.5))
pred_df$pred <- predict(fit, pred_df)

ring <- data.frame(ring_in = c(3.5, 3.75, 4)) %>% 
  mutate(ring_mm = ring_in * 25.4)

soak_times <- c(24, 48)

sel <- array(data = NA,
      dim = c(length(pred_df$pred), length(ring$ring_mm) + 1, length(soak_times)))

nsim <- 5000

for(i in 1:length(pred_df$pred)) {
  for(j in 1:length(ring$ring_mm)) {  
    for(k in 1:length(soak_times)) {  
      
      # Simulate girths based on variability (girth_se) in the girth-length
      # relationship (lognormally distributed)
      sim <- exp(rnorm(n = nsim, mean = pred_df$pred[i], sd = girth_se)) / pi
      
      # Determine which fish in the simulation are greater than the escape ring (1
      # = greater than escape ring, would be retained by pot; 0 = smaller than
      # escape ring and could escape)
      tmp <- as.numeric(sim > ring$ring_mm[j])
      escapees <- tmp[tmp == 0]
      
      # Use the probability of the escapees finding the escape ring
      # given soak time to determine fraction of potential escapees that will
      # remain in the pot due to soak time.
      prob <- soak %>% filter(hr == soak_times[k]) %>% pull(p)
      n <- round(prob * length(escapees), 0)
      escapees[(n+1):length(escapees)] <- 1
      
      # If there are no potential escapees, assign selectivity as 1
      if(n == 0) {
        sel[i,j,k] <- 1 # equivalent to sum(tmp[tmp == 1]) / nsim
      } else {
        sel[i,j,k] <- sum(tmp[tmp == 1], escapees) / nsim
      }
      
      # Soak time column
      sel[,4,k] <- soak_times[k]
    }
  }
}

sel <- rbind(as.data.frame(sel[,,1]), as.data.frame(sel[,,2]))
names(sel) <- c(levels(bio$Treatment)[2:4], "soak_time")
sel <- sel %>% mutate(length = rep(pred_df$length, length(soak_times)))
sel <- data.table::melt(data = sel, id.vars = c("soak_time", "length"), variable.name = "Treatment", value.name = "p")

# Fish below a certain size can fit through the mesh. Adjust selectivity
# obtained from including soak times accordingly
ctl <- ctl %>% filter(length %in% seq(30, 100, 0.5))

sel %>% 
  # filter(soak_time == 24) %>% 
  left_join(ctl %>% select(length, ctl_p = p), by = "length") %>% 
  mutate(p = p * ctl_p) -> sel

p <- ggplot() +
  geom_line(data = ctl, aes(x = length, y = p, lty = Treatment), colour = "grey90", size = 1) +
  geom_line(data = sel, aes(x = length, y = p, col = factor(soak_time), 
                     linetype = Treatment, group = interaction(Treatment, soak_time)), size = 0.5) +
  scale_colour_manual(values = c("grey70", "black")) +
  scale_linetype_manual(values = c(1, 2, 3, 1)) +
  labs(x = "\nLength (cm)", y = "Proportion retained\n", color = "Soak time (hr)") +
  theme(legend.position = c(0.8, 0.45),
        legend.key.width=unit(1.5,"line")) +
  guides(linetype = guide_legend(override.aes = list(size = c(0.5, 0.5, 0.5, 1),
                                                     colour = c("black", "black", "black", "grey90"))))
        # legend.key.width=unit(0.5,"cm"))
p

# Caption: Theoretical selectivity curves for escape ring treatments as a function of soak time.
ggsave(plot = p, filename = paste0("figures/theoretical_selectivity_soaktime_", YEAR, ".png"),
       dpi=300, height=3, width=6, units="in")
