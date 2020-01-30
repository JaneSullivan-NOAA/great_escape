# Exploratory data analysis
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-10-01

# Set up ----

YEAR <- 2019 # study year(s)
source("code/helper.r")
library(FSA) # for Dunn test
library(broom) # tidy

# Bio data
bio <- read_csv(paste0("data/bio_cleaned_", YEAR, ".csv")) %>% 
  mutate(Treatment = derivedFactor("Control" = Treatment == "Control",
                            "8.9 cm" = Treatment == "3.50 in",
                            "9.5 cm" = Treatment == "3.75 in",
                            "10.2 cm" = Treatment == "4.00 in",
                            .default = NA,
                            .ordered = TRUE))

# Fishery girth data (collected in order to account for grith seasonal adjustments)
fsh_grth <- read_csv(paste0("data/fsh_bio_", YEAR, ".csv")) 

# Effort data 
effort <- read_csv(paste0("data/pot_effort_", YEAR, ".csv"))

# Data sent to K. Wood for map 20191230
effort %>% 
  distinct(effort_no, start_lat, start_lon, end_lat, end_lon) %>% 
  write_csv("output/station_locations.csv")

# Total sablefish counts by depth and disposition (some pots in a set were
# dumped due to processing time)
counts <- read_csv(paste0("data/total_counts_", YEAR, ".csv")) %>% 
  mutate(Treatment = derivedFactor("Control" = treatment == "Blue",
                                   "8.9 cm" = treatment == "Purple",
                                   "9.5 cm" = treatment == "Green",
                                   "10.2 cm" = treatment == "Yellow",
                                   .default = NA,
                                   .ordered = TRUE))

# Length frequency ----

p <- bio %>% 
  filter(!is.na(Treatment) & between(length, 40, 90)) %>%
  droplevels() %>%
  ggplot(aes(x = length, colour = Treatment, linetype = Treatment)) + 
  geom_freqpoly() +
  scale_colour_manual(values = c("grey90", "grey70", "grey40", "black")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  xlim(40, 90) +
  labs(x = "Fork length (cm)", y = "Count") + 
  theme(legend.position = c(0.8, 0.65))


p
# Caption: Length frequency distribution by escape ring treatment.
ggsave(plot = p, filename = paste0("figures/size_freq_", YEAR, ".pdf"), dpi=600, height=80/1.618, width=80, units="mm")

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
  labs(x = "Fork length (cm)", y = "Girth (mm)") +
  theme(legend.position = c(0.8, 0.2))
p
ggsave(plot = p, filename = paste0("figures/girth_outliers_", YEAR, ".pdf"), 
       dpi=600, height=80/1.618, width=80, units="mm")

grth <- grth %>% filter(Outlier != "Outlier") %>% 
  select(-Outlier)

# Plot cleaned up girth measurements by treatment
p <- ggplot(grth, aes(x = length, y = girth)) +
  geom_point(shape = 20) +
  facet_wrap(~ Treatment) +
  labs(x = "Fork length (cm)", y = "Girth (mm)")
p
ggsave(plot = p, filename = paste0("figures/girth_bytreatment_", YEAR, ".pdf"), 
       dpi=600, height=180/1.618, width=180, units="mm")

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
  labs(x = "Fork length (cm)", y = "Girth (mm)") 
p
ggsave(plot = p, filename = paste0("figures/fitted_girth_bytreatment_", YEAR, ".pdf"), 
       dpi=600, height=180/1.618, width=180, units="mm")

# Girth adjustments  ----

# Combine girths from survey and fishery
comb_grth <- grth %>% 
  filter(Treatment == "Control") %>%
  select(length, girth) %>% 
  mutate(Source = "Survey (May)") %>% 
  bind_rows(fsh_grth %>% 
              select(length, girth) %>% 
              mutate(Source = "Fishery (Sep/Oct)"))

# Summary of fish girths/weights
comb_grth %>% 
  group_by(Source) %>% 
  # filter(!is.na(girth) | !is.na(weight)) %>% 
  summarise(n = n()) %>% 
  write_csv("output/girth_sample_sizes.csv")

fit_simple <- glm(log(girth) ~ log(length), 
                  family = gaussian(link = "identity"), data = comb_grth)
fit_int <- glm(log(girth) ~ log(length) + Source, 
               family = gaussian(link = "identity"), data = comb_grth)
fit_intslp <- glm(log(girth) ~ log(length) * Source, 
                  family = gaussian(link = "identity"), data = comb_grth)

AIC_simple <- AIC(fit_simple)
AIC_int <- AIC(fit_int)
AIC_intslp <- AIC(fit_intslp) # best model = different intercepts
AIC_simple - AIC_int
AIC_int - AIC_intslp

summary(fit_int)
R2_int <- 1 - (fit_int$deviance / fit_int$null.deviance) # coefficient of variation
sigma(fit_int) # residual standard deviation
# Get fitted values for girth and prediction intervals (need to exponentiate for figure)
pred <- comb_grth %>% select(Source, length)
pred <- ciTools::add_pi(pred, fit_int, alpha = 0.05, names = c("pi_lwr", "pi_upp"))

# Apply bias correction
pred <- pred %>% 
  mutate(fitted = exp(pred) * exp(0.5 * sigma(fit_int)^2),
         lower = exp(pi_lwr) * exp(0.5 * sigma(fit_int)^2),
         upper = exp(pi_upp) * exp(0.5 * sigma(fit_int)^2))

pred <- pred %>% mutate(Source = factor(Source, levels = c("Survey (May)", "Fishery (Sep/Oct)"), ordered = TRUE))
comb_grth <- comb_grth %>% mutate(Source = factor(Source, levels = c("Survey (May)", "Fishery (Sep/Oct)"), ordered = TRUE))

p1 <- ggplot() +
  geom_ribbon(data = pred, aes(x = length, ymin = lower, ymax = upper, fill = Source), 
              alpha = 0.3) +
  geom_point(data = comb_grth, aes(x = length, y = girth, colour = Source, shape = Source), size = 0.8, alpha = 0.7) +
  geom_line(data = pred, aes(x = length, y = fitted, group = Source, colour = Source, linetype = Source)) +
  scale_colour_manual(values = c("grey10", "grey60")) +
  scale_fill_manual(values = c("grey80", "grey70")) +
  labs(x = "Fork length (cm)", y = "Girth (mm)") +
  theme(legend.position = c(0.75, 0.2),
        legend.text=element_text(size = 7)) +
  annotate('text', x = 45, y = 650,
           label = as.character(expression(paste(R^{2}==0.9, ",  ", sigma==0.05))),
           parse = TRUE, size = 3, hjust = 0) +
  annotate('text', x = 45, y = 620,
           label = as.character(expression(paste("May:  ", italic(hat(G))==4.47*italic(L)^{1.03}))),
           parse = TRUE, size = 3, hjust = 0) +
  annotate('text', x = 45, y = 590,
           label = as.character(expression(paste("Sep/Oct:  ", italic(hat(G))==4.47*italic(L)^{1.07}))),
           parse = TRUE, col = "grey50", size = 3, hjust = 0)
  
p1

# Caption: A comparison of fitted values and prediction intervals for the
# regression of girth on length for data collected during the survey in May
# (black circles) and fishery in September and October (grey triangles).
ggsave(plot = p1, filename = paste0("figures/girth_bysource_", YEAR, ".pdf"), 
       dpi=600, height=180/1.618, width=180, units="mm")

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

pred_df <- data.frame(length = seq(30, 100, 0.01), 
                      Source = factor(rep("Survey (May)", length(seq(30, 100, 0.01))))) 
                                                                   
pred_df$pred <- predict(fit_int, pred_df)

# A. Baldwin measured the mesh size diameter for me. Assume 70 mm. See issue 2
# on github for documentation.
ring <- data.frame(ring_in = c(70 / 25.4, 3.5, 3.75, 4)) %>% 
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
sel <- sel %>% filter(length %in% seq(30, 100, 0.1))

# Save output, used to inform prior
write_csv(sel, paste0("output/theoretical_selectivity_", YEAR, ".csv"))

p <- ggplot(sel, aes(x = length, y = p, col = Treatment, 
                     linetype = Treatment, group = Treatment)) +
  geom_hline(yintercept = 0.5, col = "lightgrey", size = 0.4, lty = 2) +
  geom_vline(xintercept = 63, col = "lightgrey", size = 0.4, lty = 2) +
  # geom_point(aes(x = 63, y = 0.5), col = "green", size = 1) +
  geom_line() +
  scale_colour_manual(values = c("grey90", "grey70", "grey40", "black")) +
  scale_linetype_manual(values = c(1, 1, 2, 3)) +
  labs(x = "Fork length (cm)", y = "Proportion retained") +
  theme(legend.position = c(0.8, 0.3)) +
  xlim(c(35,85))
p
# Caption: Theoretical selectivity curves for control and escape ring treatments.
ggsave(plot = p, filename = paste0("figures/theoretical_selectivity_", YEAR, ".pdf"),
       dpi=600, height=80/1.618, width=80, units="mm")

# Theoretical with fishery girth ----

pred_df <- data.frame(length = seq(30, 100, 0.01), 
                      Source = factor(rep("Fishery (Sep/Oct)", length(seq(30, 100, 0.01))))) 

pred_df$pred <- predict(fit_int, pred_df)

sel_grth <- matrix(nrow = length(pred_df$pred),
              ncol = length(ring$ring_mm))
nsim <- 5000

for(i in 1:length(pred_df$pred)) {
  
  sim <- exp(rnorm(n = nsim, mean = pred_df$pred[i], sd = girth_se)) / pi
  
  for(j in 1:length(ring$ring_mm)) {
    
    sel_grth[i,j] <- length(which(sim > ring$ring_mm[j])) / nsim
  }
}

# Reformat and plot
sel_grth <- as.data.frame(sel_grth)
names(sel_grth) <- levels(bio$Treatment)
sel_grth <- sel_grth %>% mutate(length = pred_df$length)
sel_grth <- data.table::melt(data = sel_grth, id.vars = "length", variable.name = "Treatment", value.name = "p")

write_csv(sel_grth, paste0("output/theoretical_slx_fishery_", YEAR, ".csv"))

full_sel <- sel %>% 
  mutate(Source = "Survey (May)") %>% 
  bind_rows(sel_grth %>% 
              mutate(Source = "Fishery (Sep/Oct)")) %>% 
  # filter(Treatment != "Control") %>% 
  droplevels() %>% 
  mutate(Source = factor(Source, levels = c("Survey (May)", "Fishery (Sep/Oct)"), ordered = TRUE))

full_sel <- full_sel %>% filter(length %in% seq(40, 100, 0.6))

p2 <- ggplot(full_sel, aes(x = length, y = p, col = Source, 
                     linetype = Treatment, group = interaction(Source, Treatment))) +
  # geom_hline(yintercept = 0.5, col = "lightgrey", size = 0.4, lty = 2) +
  geom_vline(xintercept = 63, col = "grey85", size = 0.5, lty = 5) +
  geom_line(size = 0.7) +
  scale_colour_manual(values = c("grey10", "grey60")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  labs(x = "Fork length (cm)", y = "Proportion retained") +
  theme(legend.position = c(0.79, 0.3),
        legend.text=element_text(size = 7),
        legend.spacing.y = unit(0, "cm")) +
  # annotate("curve", x = 50, y = 0.85, xend = 62.5, yend = 1,
  #          colour = "grey70", curvature = -0.3, arrow = arrow(length = unit(1, "mm"))) +
  # annotate("text", x = 50, y = 0.8, colour = "grey60", size = 3,
  #          label = as.character(expression(paste(italic(L)[50]== "63 cm"))), parse = TRUE) +
  xlim(c(40,85))
p2
p3 <- plot_grid(p1, p2, ncol = 2, labels = c("A", "B"))
p3
ggsave(plot = p3, filename = paste0("figures/girth_regression_theoretical_slx_", YEAR, ".pdf"),
       dpi=600, height=180/1.618, width=180, units="mm")

# Values for the text showing % selected at 63 and lengths at which the
# treatments are fully selected
full_sel %>% filter(length == 63)
full_sel %>%
  filter(p > 0.999) %>% 
  group_by(Source, Treatment) %>% 
  summarize(min(length))

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
ggsave(plot = p, filename = paste0("figures/escape_prob_soaktime.pdf"),
       dpi=600, height=80/1.618, width=80, units="mm")

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
  # geom_hline(yintercept = 0.5, col = "lightgrey", size = 0.4, lty = 2) +
  # geom_vline(xintercept = 61, col = "lightgrey", size = 0.4, lty = 2) +
  geom_line(data = ctl, aes(x = length, y = p, lty = Treatment), colour = "grey90", size = 1) +
  geom_line(data = sel, aes(x = length, y = p, col = factor(soak_time), 
                     linetype = Treatment, group = interaction(Treatment, soak_time)), size = 0.5) +
  # geom_point(aes(x = 61, y = 0.5), col = "green", size = 1.5) +
  scale_colour_manual(values = c("grey70", "black")) +
  scale_linetype_manual(values = c(1, 2, 3, 1)) +
  labs(x = "Fork length (cm)", y = "Proportion retained", color = "Soak time (hr)") +
  theme(legend.position = c(0.8, 0.45),
        legend.key.width=unit(1.5,"line")) +
  guides(linetype = guide_legend(override.aes = list(size = c(0.5, 0.5, 0.5, 1),
                                                     colour = c("black", "black", "black", "grey90"))))

p

# Caption: Theoretical selectivity curves for escape ring treatments as a function of soak time.
ggsave(plot = p, filename = paste0("figures/theoretical_selectivity_soaktime_", YEAR, ".pdf"),
       dpi=600, height=180/1.618, width=180, units="mm")

# Capture efficieny ----

# Shapiro-Wilk test suggests CPUE (n salbefish per pot) is not normally
# distributed, so we did a Kruskal-Wallis test, which is essentially a
# nonparameteric one-way ANOVA
shapiro.test(counts$n_sablefish) # H0: Data are normal
kw <- tidy(kruskal.test(n_sablefish ~ Treatment, data = counts)) %>% 
  mutate(data = "All sizes combined") # H0: means of the groups are the same
kw
tst_counts <- counts %>% mutate(Treatment = factor(Treatment, ordered = FALSE))
# Use Dunn test for adhoc multiple comparisons with a Bonferroni adjusted
# p-value (P.adj)
dunn <- dunnTest(n_sablefish ~ Treatment, data = tst_counts, method = "bonferroni")
dunn_df <- as.data.frame(dunn[[2]]) %>% mutate(data = "All sizes combined")

counts %>% 
  group_by(Treatment) %>% 
  summarize(mean_cpue = mean(n_sablefish),
            median_cpue = median(n_sablefish),
            cpue_sd = sd(n_sablefish)) %>% 
  kable()

# Split up CPUE data by length. We have a smaller sample size than the combined
# Two categories: 1) < 63 cm (L50, Dressel 2009) and 2) >= 63 cm

# Total sample sizes
counts %>% 
  group_by(Treatment) %>% 
  dplyr::summarise(n = sum(n_sablefish)) %>% 
  bind_rows(counts %>% 
              mutate(Total = "Total") %>% 
              group_by(Total) %>% 
              summarize(n = sum(n_sablefish)) %>% 
              rename(Treatment = Total)) %>%
  left_join(bio %>% # Mean size overall
              group_by(Treatment) %>% 
              summarize(n_lengths = n(),
                        mean_length = mean(length),
                        se_length = sd(length)/sqrt(n_lengths)) %>% 
              bind_rows(bio %>% 
                          mutate(Total = "Total") %>% 
                          group_by(Total) %>% 
                          summarize(n_lengths = n(),
                                    mean_length = mean(length),
                                    se_length = sd(length)/sqrt(n_lengths)) %>% 
                          rename(Treatment = Total))) %>% 
  write_csv("output/sample_sizes.csv")

size_cpue <- bio %>% 
  # Remove "stuck" fish b/c they cannot be attributed to a specific pot
  filter(pot_no != 99) %>% 
  mutate(Size_category = ifelse(length < 63, 
                                "Sablefish < 63 cm", "Sablefish \u2265 63 cm")) %>%
  group_by(Treatment, effort_no, pot_no, Size_category) %>% 
  summarize(n_sablefish = n()) 

# Figure with all sizes combined, and split by L50

mu <- counts %>% 
  filter(Treatment == "Control") %>% 
  summarize(mu = median(n_sablefish)) %>% 
  mutate(Size_category = "All sizes combined") %>% 
  bind_rows(size_cpue %>% 
              filter(Treatment == "Control") %>% 
              group_by(Size_category) %>% 
              summarize(mu = median(n_sablefish)))

counts %>% 
  select(Treatment, effort_no = set, pot_no = pot, n_sablefish) %>% 
  mutate(Size_category = "All sizes combined") %>% 
  bind_rows(size_cpue) -> tmp

p <- ggplot(tmp, aes(x = Treatment, y = n_sablefish)) +
  # If the notches don't overlap this suggests the means are different
  geom_boxplot(notch = TRUE) +
  geom_hline(data = mu, aes(yintercept = mu), col = "grey", lty = 2) +
  facet_wrap(~ Size_category, ncol = 1, scales = "free_y") + 
  labs(x = NULL, y = "CPUE")
p

ggsave(filename = paste0("figures/cpue_", YEAR, ".pdf"),
       dpi=600, height=3*(80/1.618), width=80, units="mm", device = cairo_pdf)

# All significant
tst <- size_cpue %>% 
  ungroup() %>% 
  filter(Size_category == "Sablefish < 63 cm") %>% 
  mutate(Treatment = factor(Treatment, ordered = FALSE))
tmp <- tidy(kruskal.test(n_sablefish ~ Treatment, data = tst)) %>% 
  mutate(data = "Sablefish < 63 cm")# H0: means of the groups are the same
kw <- kw %>% bind_rows(tmp)
dunn <- dunnTest(n_sablefish ~ Treatment, data = tst, method = "bonferroni") # Bonferroni
dunn
dunn_df <- dunn_df %>% 
  bind_rows(as.data.frame(dunn[[2]]) %>% 
              mutate(data = "Sablefish < 63 cm"))

tst %>% 
  group_by(Treatment) %>% 
  summarize(median = median(n_sablefish),
            mean = mean(n_sablefish))

tst <- size_cpue %>% ungroup() %>% 
  filter(Size_category == "Sablefish \u2265 63 cm") %>% 
  mutate(Treatment = factor(Treatment, ordered = FALSE))
tmp <- tidy(kruskal.test(n_sablefish ~ Treatment, data = tst)) %>% 
  mutate(data = "Sablefish >= 63 cm")# H0: means of the groups are the same
kw <- kw %>% bind_rows(tmp)
dunn <- dunnTest(n_sablefish ~ Treatment, data = tst, method = "bonferroni") # Bonferroni
dunn
dunn_df <- dunn_df %>% 
  bind_rows(as.data.frame(dunn[[2]]) %>% 
              mutate(data = "Sablefish >= 63 cm"))

size_cpue %>% ungroup() %>% 
  filter(Size_category == "Sablefish \u2265 63 cm") %>% 
  group_by(Treatment) %>% 
  summarize(median = median(n_sablefish),
            mean = mean(n_sablefish))

write_csv(kw, "output/kw_tests.csv")
write_csv(dunn_df, "output/dunn_tests.csv")
