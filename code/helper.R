# Libraries, ggplot themes, and user-defined fxns.
# Contact: jane.sullivan@noaa.gov
# Last edited: July 2020

options(scipen = 999) # turn off scientific notation

# Libraries ----

if(!require("mosaic"))   install.packages("mosaic") # derivedFactor, derivedVariable. Masks over a lot of fxns, but generally only improves their utility
if(!require("tidyverse"))   install.packages("tidyverse") # dplyr, ggplot, etc.
if(!require("lubridate"))   install.packages("lubridate") # dates functions like yday, dmy/mdy
if(!require("data.table"))   install.packages("data.table") # dcast, foverlaps
# if(!require("ROracle"))   install.packages("ROracle") # database access through R
if(!require("knitr"))   install.packages("knitr") # r markdown
if(!require("cowplot"))   install.packages("cowplot") # plot_grid and so much else
# if(!require("captioner"))   install.packages("captioner") #numbering, ordering, & creating captions for tables and figures
if(!require("ggpubr"))   install.packages("ggpubr") # QQ plots
if(!require("ciTools"))   install.packages("ciTools") # GLM prediction intervals
if(!require("FSA"))   install.packages("FSA") # Dunn tests
if(!require("broom"))   install.packages("broom") # tidy GLM results
if(!require("boot"))   install.packages("boot") # bootstrapping
if(!require("xtable"))   install.packages("xtable") # generate latex tables

# Figure theme ----

windowsFonts(Times=windowsFont("Helvetica"))

theme_sleek <- function(base_size = 13, base_family = "Helvetica") {
  half_line <- base_size/2
  theme_light(base_size = 13, base_family = "Helvetica") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      #axis.text = element_text(colour = "grey30"),
      #axis.title = element_text(colour = "grey30"),
      # legend.title = element_text(size = rel(0.9)),
      panel.border = element_rect(fill = NA),#, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      # legend.text = element_text(size = 11),
      legend.key.height = unit(0.6, 'line'),
      # legend.spacing.y = unit(0, "cm"),
      # legend.text = element_text(size = rel(0.7)),#, colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)#,
      #plot.title = element_text(colour = "grey30"),#, size = rel(1)
      #plot.subtitle = element_text(colour = "grey30")#, size = rel(.85)
    )
}
theme_set(theme_sleek())

# User-defined functions ----

# Format ggplot figures with ticked axes (especially good for marking year and
# age) 

# Depends on dplyr
tickr <- function(
  data, # dataframe
  var, # column of interest
  to # break point definition 
){
  
  VAR <- enquo(var) # makes VAR a dynamic variable
  
  data %>% 
    distinct(!!VAR) %>%
    ungroup(!!VAR) %>% 
    mutate(labels = ifelse(!!VAR %in% seq(to * round(min(!!VAR) / to), max(!!VAR), to),
                           !!VAR, "")) %>%
    select(breaks = UQ(VAR), labels)
}

#fig_nums, tbl_nums, and appendix_nums fxns created fromm captioner() fxn in
#'captioner' library, which I've tweaked below, changed separator from a colon
#to period) - these fxns are used for autonumbering figs and tables in text and
#creating captions.

# captioner <- function (prefix = "Figure", auto_space = TRUE, levels = 1, 
                       # type = NULL, infix = ".") {
#   check_class(prefix, "character")
#   check_class(auto_space, "logical")
#   check_class(levels, "numeric")
#   check_class(infix, "character")
#   if (is.null(type)) {
#     type <- c(rep("n", times = levels))
#   }
#   else if (length(type) < levels) {
#     type[(length(type) + 1):levels] <- "n"
#   }
#   else if (length(type) > levels) {
#     type <- type[1:levels]
#   }
#   if (!all(type %in% c("n", "c", "C"))) {
#     stop("Invalid 'type' value used.  Expecting 'n', 'c', or 'C'.")
#   }
#   if (auto_space) {
#     prefix <- paste(prefix, " ")
#   }
#   force(levels)
#   force(prefix)
#   force(infix)
#   OBJECTS <- list(name = NULL, caption = NULL, number = list(list()))
#   OBJECTS$number[[1]][which(type == "n")] <- 1
#   OBJECTS$number[[1]][which(type == "c")] <- "a"
#   OBJECTS$number[[1]][which(type == "C")] <- "A"
#   function(name, caption = "", display = "full", level = FALSE, 
#            cite = FALSE, num = FALSE) {
#     if (level > levels) {
#       stop("Level too large.")
#     }
#     objects <- OBJECTS
#     if (any(objects$name == name)) {
#       obj_ind <- match(name, objects$name)
#       if (objects$caption[obj_ind] == "") {
#         objects$caption[obj_ind] <- caption
#       }
#       else {
#         caption <- objects$caption[obj_ind]
#       }
#     }
#     else {
#       obj_ind <- length(objects$name) + 1
#       if (length(objects$number) == length(objects$name)) {
#         if (level) {
#           objects$number[[obj_ind]] <- increment(objects$number[[obj_ind - 
#                                                                    1]], level)
#         }
#         else {
#           objects$number[[obj_ind]] <- increment(objects$number[[obj_ind - 
#                                                                    1]], levels)
#         }
#       }
#       objects$name[obj_ind] <- name
#       objects$caption[obj_ind] <- caption
#     }
#     assign("OBJECTS", objects, envir = parent.env(environment()))
#     obj_num <- paste(objects$number[[obj_ind]], collapse = infix)
#     if (cite) {
#       .Deprecated(new = "display", old = "cite")
#       return(paste0(prefix, obj_num))
#     }
#     if (num) {
#       .Deprecated(new = "display", old = "num")
#       return(obj_num)
#     }
#     if (display == FALSE) {
#       return(invisible())
#     }
#     #FLAG: Jane changed ": " to ". "
#     else if (display == "full" || display == "f") {
#       return(paste0(prefix, obj_num, ". ", caption))
#     }
#     else if (display == "cite" || display == "c") {
#       return(paste0(prefix, obj_num))
#     }
#     else if (display == "num" || display == "n") {
#       return(obj_num)
#     }
#     else {
#       warning("Invalid display mode used.  Caption was still saved.")
#       return(invisible())
#     }
#   }
# }

# fig <- captioner(prefix = "Figure")
# tbl <- captioner(prefix = "Table")
# appendix_tbl <- captioner(prefix = "Table") #Numbers tables in the appendix
# appendix_fig <- captioner(prefix = "Figure") #Numbers figures in the appendix

#
# Calculate marginal AIC for a fitted model
# From: https://github.com/kaskr/TMB_contrib_R/blob/master/TMBhelper/R/TMBAIC.R (J. Thorson?)

TMBAIC = function(opt, # output from nlminb or optim
                p = 2, # penalty on additional fixed effects (default=2, for AIC)
                n = Inf) { # sample size, for use in AICc calculation (default=Inf, for which AICc=AIC)
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
  return( Return )
}

# Plot MLE results
plot_mle <- function(delta = c("estimated", "fixed"), # Model version with estimated or fixed delta
                     selprior = c("beta", "normal"), # prior on selectivity at 0, 1
                     save = TRUE) {
  
  phi <- as.data.frame(obj$report()$fit_phi)
  names(phi) <- paste0(unique(df$Treatment))
  phi <- phi %>% 
    mutate(length_bin = len) %>% 
    pivot_longer(-length_bin, names_to = "Treatment", values_to = "phi")
  # phi <- data.table::melt(data = phi, id.vars = "length_bin", variable.name = "Treatment", value.name = "phi") 
  
  com %>% 
    left_join(phi, by = c("Treatment", "length_bin")) %>% 
    mutate(resid = combined_p - phi,
           Treatment = factor(Treatment, levels = c("8.9 cm", "9.5 cm", "10.2 cm"), 
                              ordered = TRUE)) -> tmp
  
  p1 <- ggplot(tmp) +
    geom_point(aes(x = length_bin, y = combined_p,
                   group =  Treatment, col = Treatment)) +
    geom_line(aes(x = length_bin, y = phi, 
                  group = Treatment, col = Treatment)) +
    labs(x = "Length", y = "Proportion in treatment pots") 
  
  
  p2 <- ggplot(tmp, aes(x = length_bin, y = resid, col = Treatment)) +
    geom_hline(yintercept = 0, col = "grey") +
    geom_segment(aes(x = length_bin, xend = length_bin, y = 0, 
                     yend = resid),
                 size = 0.2, col = "grey") +
    geom_point() +
    facet_wrap(~Treatment) +
    theme(legend.position = "none") +
    labs(x = "Length", y = "Residuals")
  
  
  # Plot selectivity
  slx <- as.data.frame(obj$report()$full_slx)
  names(slx) <- paste0(unique(df$Treatment))
  slx <- slx %>% 
    mutate(length_bin = fit_len) %>% 
    pivot_longer(-length_bin, names_to = "Treatment", values_to = "slx") %>% 
    mutate(Treatment = derivedFactor("Control" = Treatment == "Control",
                                     "8.9 cm" = Treatment == "8.9 cm",
                                     "9.5 cm" = Treatment == "9.5 cm",
                                     "10.2 cm" = Treatment == "10.2 cm",
                                     .default = NA,
                                     .ordered = TRUE)) 
  
  p3 <- ggplot() +
    ylim(c(0,1)) +
    geom_line(data = slx, aes(x = length_bin, y = slx, 
                              group = Treatment, col = Treatment)) +
    labs(x = "Length", y = "Proportion selected") +
    theme(legend.position = "none")
  
  p4 <- data.frame(nu = best[8:(length(best))]) %>% 
    ggplot(aes(x = nu)) +
    geom_histogram(bins = 15, fill = "white", col = "black") +
    labs(x = "Set-level random effects", y = NULL)
  
  pa <- plot_grid(p1, p2, ncol = 1,rel_heights = c(1/3, 2/3))
  pb <- plot_grid(p3, p4, ncol = 1, rel_heights = c(.6, .4))
  pab <- plot_grid(pa, pb, rel_widths = c(.6, .4))
  
  # now add the title
  title <- ggdraw() + 
    draw_label(paste0("MLE results (delta = ", delta,
                      "; prior on selectivity 0 and 1 = ", selprior, ")"),
               fontface = 'bold', x = 0, hjust = 0) +
   
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
       theme(plot.margin = margin(0, 0, 0, 7))
  finalp <- plot_grid(title, pab, ncol = 1, rel_heights = c(0.1, 1))
  
  # save fig
  if(save == TRUE) {
    ggsave(plot = finalp,
           filename = paste0("../figures/mle_fit_", delta, "_", selprior, ".pdf"),
           width = 180, height = 180/1.618, units = "mm")
  }
  
  print(finalp)
}
