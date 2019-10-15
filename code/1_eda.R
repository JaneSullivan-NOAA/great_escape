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

# Effort data 
effort <- read_csv(paste0("data/pot_effort_", YEAR, ".csv"))

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

# Girth adjustments  ----

# Lay groundwork for girth adjustments - currently not used

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
pred <- ciTools::add_pi(pred, fit_complex, alpha = 0.5, names = c("pi_lwr", "pi_upp"))

p <- ggplot() +
  geom_point(data = grth, aes(x = length, y = girth), shape = 20) +
  geom_ribbon(data = pred, aes(x = length, ymin = exp(pi_lwr), ymax = exp(pi_upp)), 
              alpha = 0.6, fill = "grey70") +
  geom_line(data = pred, aes(x = length, y = exp(pred), group = Treatment)) +
  facet_wrap(~ Treatment) +
  labs(x = "\nLength (cm)", y = "Girth (mm)\n") 

ggsave(plot = p, filename = paste0("figures/fitted_girth_bytreatment_", YEAR, ".png"), 
       dpi=300, height=5, width=6, units="in")

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

pred_df <- data.frame(length = seq(30, 100, 1))
pred_df$pred <- predict(fit, pred_df)

# Assume mesh size diam of pots = 2.5 in until I can get better info from A.
# Baldwin (email sent 2019-10-04)
ring <- data.frame(ring_in = c(2.5, 3.5, 3.75, 4)) %>% 
  mutate(ring_mm = ring_in * 25.4)

# Assume girths are lognormally distributed. Simulate girth distribution at 1 cm
# increments, divide by pi to get approximate fish girth. Determine proportion
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

p <- ggplot(sel, aes(x = length, y = p, col = Treatment, 
                     linetype = Treatment, group = Treatment, size = Treatment)) +
  geom_line() +
  scale_colour_manual(values = c("grey90", "grey70", "grey40", "black")) +
  scale_size_manual(values = c(1.5, 0.7, 0.7, 0.7)) +
  scale_linetype_manual(values = c(1, 1, 2, 3)) +
  labs(x = "\nLength (cm)", y = "Proportion retained\n") +
  theme(legend.position = c(0.8, 0.3))

# Caption: Theoretical selectivity curves for control and escape ring treatments.
ggsave(plot = p, filename = paste0("figures/theoretical_selectivity_", YEAR, ".png"),
       dpi=300, height=3, width=6, units="in")

# Save output, 90% selectivity used to inform prior
write_csv(sel, paste0("output/theoretical_selectivity_", YEAR, ".csv"))
