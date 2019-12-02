# Libraries, ggplot themes, and user-defined fxns.
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-10-01

options(scipen = 999) # turn off scientific notation

# Libraries ----

if(!require("mosaic"))   install.packages("mosaic") # derivedFactor, derivedVariable. Masks over a lot of fxns, but generally only improves their utility
# if(!require("tidyverse"))   install.packages("tidyverse") # dplyr, ggplot, etc.
if(!require("dplyr"))   install.packages("dplyr") # dplyr, ggplot, etc.
if(!require("ggplot2"))   install.packages("ggplot2") # dplyr, ggplot, etc.
if(!require("readr"))   install.packages("readr") # read_csv
# if(!require("tidyr"))   install.packages("tidyr") 
if(!require("lubridate"))   install.packages("lubridate") # dates functions like yday, dmy/mdy
# if(!require("mgcv"))   install.packages("mgcv") # gams
# if(!require("gridExtra"))   install.packages("gridExtra") # multipanneled plots
if(!require("data.table"))   install.packages("data.table") # dcast, foverlaps
if(!require("ROracle"))   install.packages("ROracle") # database access through R
# if(!require("broom"))   install.packages("broom") # tidying regression model output
# if(!require("padr"))   install.packages("padr") # fills in missing values in a time series
# if(!require("tidyr"))   install.packages("tidyr") # reshaping data
if(!require("knitr"))   install.packages("knitr") # r markdown
# if(!require("forcats"))   install.packages("forcats") # releveling factors
# if(!require("cowplot"))   install.packages("cowplot") # plot_grid and so much else
# if(!require("ggridges"))   install.packages("ggridges") # length comps
if(!require("captioner"))   install.packages("captioner") #numbering, ordering, & creating captions for tables and figures
if(!require("ggpubr"))   install.packages("ggpubr") # QQ plots
if(!require("ciTools"))   install.packages("ciTools") # GLM prediction intervals

# Figure theme ----

windowsFonts(Times=windowsFont("Helvetica"))

theme_sleek <- function(base_size = 10, base_family = "Helvetica") {
  half_line <- base_size/2
  theme_light(base_size = 10, base_family = "Helvetica") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      #axis.text = element_text(colour = "grey30"),
      #axis.title = element_text(colour = "grey30"),
      #legend.title = element_text(colour = "grey30"),#, size = rel(0.9)
      panel.border = element_rect(fill = NA),#, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      #legend.text = element_text(size = rel(0.7)),#, colour = "grey30"),
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

captioner <- function (prefix = "Figure", auto_space = TRUE, levels = 1, 
                       type = NULL, infix = ".") {
  check_class(prefix, "character")
  check_class(auto_space, "logical")
  check_class(levels, "numeric")
  check_class(infix, "character")
  if (is.null(type)) {
    type <- c(rep("n", times = levels))
  }
  else if (length(type) < levels) {
    type[(length(type) + 1):levels] <- "n"
  }
  else if (length(type) > levels) {
    type <- type[1:levels]
  }
  if (!all(type %in% c("n", "c", "C"))) {
    stop("Invalid 'type' value used.  Expecting 'n', 'c', or 'C'.")
  }
  if (auto_space) {
    prefix <- paste(prefix, " ")
  }
  force(levels)
  force(prefix)
  force(infix)
  OBJECTS <- list(name = NULL, caption = NULL, number = list(list()))
  OBJECTS$number[[1]][which(type == "n")] <- 1
  OBJECTS$number[[1]][which(type == "c")] <- "a"
  OBJECTS$number[[1]][which(type == "C")] <- "A"
  function(name, caption = "", display = "full", level = FALSE, 
           cite = FALSE, num = FALSE) {
    if (level > levels) {
      stop("Level too large.")
    }
    objects <- OBJECTS
    if (any(objects$name == name)) {
      obj_ind <- match(name, objects$name)
      if (objects$caption[obj_ind] == "") {
        objects$caption[obj_ind] <- caption
      }
      else {
        caption <- objects$caption[obj_ind]
      }
    }
    else {
      obj_ind <- length(objects$name) + 1
      if (length(objects$number) == length(objects$name)) {
        if (level) {
          objects$number[[obj_ind]] <- increment(objects$number[[obj_ind - 
                                                                   1]], level)
        }
        else {
          objects$number[[obj_ind]] <- increment(objects$number[[obj_ind - 
                                                                   1]], levels)
        }
      }
      objects$name[obj_ind] <- name
      objects$caption[obj_ind] <- caption
    }
    assign("OBJECTS", objects, envir = parent.env(environment()))
    obj_num <- paste(objects$number[[obj_ind]], collapse = infix)
    if (cite) {
      .Deprecated(new = "display", old = "cite")
      return(paste0(prefix, obj_num))
    }
    if (num) {
      .Deprecated(new = "display", old = "num")
      return(obj_num)
    }
    if (display == FALSE) {
      return(invisible())
    }
    #FLAG: Jane changed ": " to ". "
    else if (display == "full" || display == "f") {
      return(paste0(prefix, obj_num, ". ", caption))
    }
    else if (display == "cite" || display == "c") {
      return(paste0(prefix, obj_num))
    }
    else if (display == "num" || display == "n") {
      return(obj_num)
    }
    else {
      warning("Invalid display mode used.  Caption was still saved.")
      return(invisible())
    }
  }
}

fig <- captioner(prefix = "Figure")
tbl <- captioner(prefix = "Table")
appendix_tbl <- captioner(prefix = "Table") #Numbers tables in the appendix
appendix_fig <- captioner(prefix = "Figure") #Numbers figures in the appendix

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
