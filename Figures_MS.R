################################################################################
# Script to reproduce Figures in the manuscript
# (note that these were edited in the paper to add silhouettes and text)
# Authors: Guillermo García-Gómez (guillegar.gz@gmail.com) & Matthew Spencer
# Date: 08062023
# Operating System: Windows 10 Pro 21H2
# ------------------------------------------------------------------------------
# Cite as: ?
# ------------------------------------------------------------------------------
rm(list=ls())# clear the work environment
# ------------------------------------------------------------------------------
# It needs to be set to Project directory
getwd()# to check
#
# ------------------------------------------------------------------------------
# Packages ####
library(brms)
library(tidyverse)
library(tidybayes)
library(plyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(lmerTest)
library(dplyr)
library(broom)
library(scales)

#### Load data ####
all_temp <- read.csv("table_S1.csv") # Table S1
all_act <- read.csv("table_S2.csv") # Table S2

# Filter data sets by group (water- or air-breather):

# Temperature-increased L:
w_temp <- all_temp %>%
  filter(group == "water-breather")

a_temp <- all_temp %>%
  filter(group == "air-breather")

# Activity-increased L:
w_act <- all_act %>%
  filter(group == "water-breather")

a_act <- all_act %>%
  filter(group == "air-breather")

# ------------------------------------------------------------------------------

# Figure 1 ####

# Create a function for the variation in the scaling slope b with temperature and activity #

# This is a conceptual schematic based on equations (2) and (3) given in Glazier (2010), 
# to illustrate the typical predictions of the Metabolic Level Boundaries Hypothesis

simulatedata <- function(n,
                         b_m, b_SA, 
                         LT0, Q10,
                         Trange, T0, T1, T2,
                         b_a, AS, act_L) {
  
  # Generate a log sequence of temperature values:
  
  # Temperature range
  temp <- seq(from = Trange[1], to = Trange[2], length.out = n)
  
  # Model metabolic level L, whose temperature dependence is determined by Q10
  # i.e., L increases exponentially with temperature
  L <- (LT0 * Q10 ^ ((temp - T0) / 10))
  
  # Modelling slope b:
  
  # Temperature effects #
  p <- seq(from = 1, to = 0, length.out = n) # relative influence of V-related processes in b (0: no influence, 1, total influence)
  
  # slope b at rest (br) as influenced by SA-related processes (approaching b_SA), 
  # and V-related processes (approaching b_m). The relative influence of these two parameters
  # is given by p:
  br <- (1 - p) * (b_SA) + p * (b_m) 
  
  ## Activity effects ##
  # generate a sequence of activity levels:
  act_L <- seq(from = Act_level[1], to = Act_level[2], length.out = n)
  
  # calculate three starting temperatures (cold, reference, hot):
  Lr_T1 = LT0 * Q10 ^ ((T1 - T0) / 10) # resting L at cold temperature (T1)
  Lr_T0 = LT0 * Q10 ^ ((T0 - T0) / 10) # resting L at reference temperature (T0)
  Lr_T2 = LT0 * Q10 ^ ((T2 - T0) / 10) # resting L at hot temperature (T2)
  
  # Assuming that the contribution of V-related processes (p) to L decreases with temperature
  # (p decreases from cold (T1) to hot temperature (T2)):
  p_T1 <- 0.9 
  p_T0 <- 0.5
  p_T2 <- 0.1
  
  # Calculate slope b at rest (br) at different temperatures:
  br_T1 <- (1 - p_T1) * (b_SA) + p_T1 * (b_m) # slope b at cold temperature (high influence of V-related processes)
  br_T0 <- (1 - p_T0) * (b_SA) + p_T0 * (b_m) # slope b at reference temperature (similar influence of V- and SA-related processes)
  br_T2 <- (1 - p_T2) * (b_SA) + p_T2 * (b_m) # slope b at hot temperature (high influence of SA-related processes)
  
  # Model activity #
  # Muscular activity leads to an exponential increase in L from rest to maximum activity 
  # (i.e., as L approximates max. aerobic capacity AS):
  L_T1 = Lr_T1 * AS ^ (act_L / 1) # total L at cold temperature (T1) including activity
  L_T0 = Lr_T0 * AS ^ (act_L / 1) # total L at reference temperature (T0) including activity
  L_T2 = Lr_T2 * AS ^ (act_L / 1) # total L at hot temperature (T2) including activity
  
  # Model the relative influence of (V-related) activity in the slope b: 
  # 0: no influence, 1: V-related processes due to activity totally dictate the slope b
  z1 <- seq(from = 0, to = 1, length.out = n)
  z2 <- seq(from = 0, to = 1, length.out = n)
  z3 <- seq(from = 0, to = 1, length.out = n)
  
  # Calculate the slope b during activity (ba) at different temperatures:
  ba_T1 <- (1 - z1) * (br_T1) + z1 * (b_a) # change in b as activity increases at cold temperature T1
  ba_T0 <- (1 - z2) * (br_T0) + z2 * (b_a) # change in b as activity increases at reference temperature T0
  ba_T2 <- (1 - z3) * (br_T2) + z3 * (b_a) # change in b as activity increases at hot temperature T2
  
  return(data.frame(L = L,
                    p = p, br = br, temp = temp,
                    act_L, p_T1 = p_T1, p_T0 = p_T0, p_T2 = p_T2,
                    br_T1 = br_T1, br_T0 = br_T0, br_T2 = br_T2,
                    z1 = z1, z2 = z2, z3 = z3,
                    ba_T1 = ba_T1, ba_T0 = ba_T0, ba_T2 = ba_T2,
                    Lr_T1 = Lr_T1, Lr_T0 = Lr_T0, Lr_T2 = Lr_T2,
                    L_T1 = L_T1, L_T0 = L_T0, L_T2 = L_T2))
}


{ # Generate values for the function #
  
  # Number of temperature values:
  n <- 1e2
  
  # Temperature effects:
  
  # Temperature range (in degrees Celsius):
  Trange <- c(0, 30)
  
  # Reference temperature (T0, in degrees Celsius):
  T0 = 15
  
  # Cold temperature (T1, in degrees Celsius):
  T1 = T0 - 12
  
  # Hot temperature (T2, in degrees Celsius):
  T2 = T0 + 12
  
  # Temperature dependence of metabolic level (Q10):
  Q10 = 2.5
  
  # Metabolic level (L) at reference temperature (T0):
  LT0 = 0.1

  # Mass-scaling (b) of V-related (i.e., maintenance) and SA-related physiological processes:
  b_m = 1 # V-related processes
  b_SA = 2/3 # SA-related processes
  
  # Activity effects:
  
  # Number of activity values:
  n2 = 1e2
  
  # Activity level ranges from 0 (no activity) to 1 (maximum activity):
  Act_level = c(0, 1)
  
  # Aerobic scope (AS):
  AS = 5 # i.e., 5 times the L at rest (this is just an illustrative example)
  
  # Mass-scaling (b) of V-related physiological processes due to muscular work:
  b_a = 1
  
}

# Simulate data using function and given values:
simdata <- simulatedata(n = n, 
                        Trange = Trange,
                        LT0 = LT0, 
                        b_SA = b_SA, 
                        b_m = b_m, 
                        b_a = b_a,  
                        Q10 = Q10, 
                        T0 = T0, T1 = T1, T2 = T2, 
                        AS = AS, act_L = act_L)

# Generate dataframe with simulated data:
b_L_total <- data.frame(cbind(temp = c(T0, T1, T2),
                              br_temp = c(unique(simdata$br_T0), unique(simdata$br_T1), unique(simdata$br_T2)),
                              Lr_temp = c(unique(simdata$Lr_T0), unique(simdata$Lr_T1), unique(simdata$Lr_T2)),
                              b_act = c(max(simdata$ba_T0),max(simdata$ba_T1),max(simdata$ba_T2)),
                              L_act = c(max(simdata$L_T0), max(simdata$L_T1), max(simdata$L_T2))))


# Generate Figure 1:

#  Change in metabolic scaling slope (b) vs metabolic level (L) 
#  under different temperature and activity level:

(Fig_1 <- ggplot(simdata, aes(x = log(L), y = br)) +  #
    geom_line(lwd = 2, lty = 1, aes(col = temp))+
    geom_line(aes(x = log(L_T0), y = ba_T0), lwd=1.25, lty=2, col = "grey70")+ 
    geom_line(aes(x = log(L_T1), y = ba_T1), lwd=1.25, lty=2, col = "grey70")+
    geom_line(aes(x = log(L_T2), y = ba_T2), lwd=1.25, lty=2, col = "grey70")+
    geom_point(data = b_L_total, aes(x = log(Lr_temp), y = br_temp, col = temp), alpha = 1, size = 6, pch = 16) + #
    geom_point(data = b_L_total, aes(x = log(L_act), y = b_act, col = temp), alpha = 1, size = 6, pch = 16) + #
    geom_point(data = b_L_total, aes(x = log(Lr_temp), y = br_temp), alpha = 1, size = 6, pch = 1) + #
    geom_point(data = b_L_total, aes(x = log(L_act), y = b_act), alpha=1, size = 6, pch = 1) + #
    scale_colour_gradient2(low = "blue", mid ="lightblue", high = "red", midpoint = 20) + # Temperature
    geom_hline(yintercept = c(2/3,1),linetype = "dashed") +
    scale_x_continuous(bquote(""~log~"Metabolic level ("~italic(L)~")"), 
                       limits=c(min(log(simdata$L))-0.25, max(log(simdata$L_T2))+0.3)) +
    scale_y_continuous("Metabolic scaling slope ( "~italic(b)~")", breaks=seq(0.4,1.6,0.1), limits=c(0.57,1.1))+
    labs(expression("Temperature ("*degree*C*")")) +
    theme_classic() +
    theme(legend.position = "none", # legend is not shown in last version of the figure
          axis.title = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.ticks.length = unit(.25, "cm"),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 13),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

# ------------------------------------------------------------------------------

# Figure 2 ####

# Load brms models:
m_all.temp_complex.info_saved <- readRDS("../Model_outputs/m_all.temp.rds") # temperature effect
m_all.act_complex.info_saved <- readRDS("../Model_outputs/m_all.act.rds") # activity effect

# get posterior distributions of brms models:
post_temp <- posterior_samples(m_all.temp_complex.info_saved) # temperature effect
post_act <- posterior_samples(m_all.act_complex.info_saved) # activity effect

### Temperature ###

# get mean estimated effects and 95% credible intervals for each group (water-, air-breathers):

post_temp %>%
  transmute(gamma_water    = b_log10_L + `b_log10_L:groupwaterMbreather`,
            gamma_notwater = b_log10_L) %>%
  gather(key, value) %>%
  group_by(key) %>%
  dplyr::summarise(mean = mean(value), lower_95ci = quantile(value, probs = c(0.025)),
                   upper_95ci = quantile(value, probs = c(0.975)))


# Plot water-breathing species data and mean estimate for log10 temperature-increased L #

# Mean estimate of log10-L effect and posterior distribution:

(wt_logL_estimate <- 
    post_temp %>%
    transmute(gamma_water    = b_log10_L + `b_log10_L:groupwaterMbreather`) %>% # get effect of Temperature-increased L in water-breathers from the interaction term
    gather(key, value) %>%
    
    # plot effect:
    ggplot(aes(x = value, group = key, color = key, fill = key)) +
    scale_x_continuous(expression(beta[italic(L)]), expand = c(0, 0), limits = c(-0.44, 0.447)) +
    scale_y_continuous(NULL, breaks = NULL) +
    stat_halfeye(aes(fill = stat(cut_cdf_qi(cdf,.width = c(1, .95)))), 
                 col = "black", .width = 0.95, point_size = 3, slab_color = "cornflowerblue") +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
    theme_bw() + 
    theme(legend.position = "none", plot.background = element_rect(fill ="transparent", color = NA)))


# Plot data of temperature-increased L in water-breathers:

(plot_wtemp <- ggplot(w_temp, aes(x = log10(L), y = b, group = experiment)) +  # Set up canvas with outcome variable on y-axis
    geom_line(color="gray60")+
    geom_point(aes(color = temp), size = 3, stroke = 1.5, alpha = 0.5) + # Plot the actual points
    scale_color_gradient2(low = "blue", mid ="lightblue", high = "red", midpoint = 17, 
                          limits = c(min(all_temp$temp), max(all_temp$temp))) + # Temperature code: colour mapped here
    geom_hline(yintercept = c(2/3,1),linetype = "dashed") +
    scale_x_continuous(bquote(''~log[10]~'Metabolic level' ~italic(L)~ '('*'mg'~O[2]~ g^-1~h^-1*')'), 
                       breaks=seq(-2.5, 0.5, 1), limits = c(-2.6,0.8))+
    scale_y_continuous("Slope"  ~italic(b)~ "", breaks = seq(0.2,1.6,0.2), limits = c(0.20,1.6))+
    labs(col = expression("T ("*degree*C*")")) +
    theme_classic() +
    theme(#legend.position = "none",
      axis.title = element_text(size = 15),
      axis.ticks.length = unit(.25, "cm"),
      axis.text = element_text(size = 13),
      axis.line = element_line(size = 0.75),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()))

# altogether:

(Fig_2A <- plot_wtemp + theme(text = element_text(size = 13),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
    annotation_custom(ggplotGrob(wt_logL_estimate), xmin = -0.55, xmax = 0.95,  
                      ymin = 1.20, 
                      ymax = 1.67))


# Plot air-breathing species data and mean estimate for log10 temperature-increased L #

# Mean estimate of log10-L effect and posterior distribution:

(at_logL_estimate <- 
    post_temp %>%
    transmute(gamma_nowater    = b_log10_L) %>% # get effect of Temperature-increased L in air-breathers from the interaction term
    gather(key, value) %>%
    
    # plot effect:
    ggplot(aes(x = value, group = key, color = key, fill = key)) +
    scale_x_continuous(expression(beta[italic(L)]), expand = c(0, 0), limits = c(-0.44, 0.447)) +
    scale_y_continuous(NULL, breaks = NULL) +
    stat_halfeye(aes(fill = stat(cut_cdf_qi(cdf,.width = c(1, .95)))), 
                 col = "black", .width = 0.95, point_size = 3, slab_color = "chartreuse3") +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_brewer(palette = "Greens", direction = -1, na.translate = FALSE) +
    theme_bw() + 
    theme(legend.position = "none", plot.background = element_rect(fill= "transparent", color = NA)))


# Plot data of temperature-increased L in air-breathers:

(plot_atemp <- ggplot(a_temp, aes(x = log10(L), y = b, group = experiment)) +  # Set up canvas with outcome variable on y-axis
    geom_line(color="gray60")+
    geom_point(aes(color = temp), size=3, stroke=1.5,alpha=0.5) + # Plot the actual points
    scale_color_gradient2(low = "blue", mid="lightblue", high = "red", midpoint=17, 
                          limits = c(min(all_temp$temp), max(all_temp$temp))) + # Temperature code: colour mapped here
    geom_hline(yintercept = c(2/3,1),linetype = "dashed") +
    scale_x_continuous(bquote(''~log[10]~'Metabolic level' ~italic(L)~ '('*'mg'~O[2]~ g^-1~h^-1*')'), 
                       breaks=seq(-2.5, 0.5, 1), limits = c(-2.6,0.8))+
    scale_y_continuous("Slope"  ~italic(b)~ "", breaks = seq(0.2,1.6,0.2), limits = c(0.20,1.6))+
    labs(col = expression("T ("*degree*C*")")) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          axis.ticks.length = unit(.25, "cm"),
          axis.text = element_text(size = 13),
          axis.line = element_line(size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

# altogether:

(Fig_2B <- plot_atemp + theme(text = element_text(size=13),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
    annotation_custom(ggplotGrob(at_logL_estimate), xmin = -0.55, xmax = 0.95,  
                      ymin = 1.20, 
                      ymax = 1.67))


### Activity level ###

# get mean estimated effects and 95% credible intervals for each group (water-, air-breathers):

post_act %>%
  transmute(gamma_water    = b_log10_L + `b_log10_L:groupwaterMbreather`,
            gamma_notwater = b_log10_L) %>%
  gather(key, value) %>%
  group_by(key) %>%
  dplyr::summarise(mean = mean(value), lower_95ci = quantile(value, probs = c(0.025)),
                   upper_95ci = quantile(value, probs = c(0.975)))

# Plot water-breathing species data and mean estimate for log10 activity-increased L #

# Mean estimate of log10-L effect and posterior distribution:

(wa_logL_estimate <- 
    post_act %>%
    transmute(gamma_water    = b_log10_L + `b_log10_L:groupwaterMbreather`) %>%
    gather(key, value) %>%
    
    ggplot(aes(x = value, group = key, color = key, fill = key)) +
    scale_x_continuous(expression(beta[italic(L)]), expand = c(0, 0), limits = c(-0.44, 0.447)) +
    scale_y_continuous(NULL, breaks = NULL) +
    stat_halfeye(aes(fill = stat(cut_cdf_qi(cdf,.width = c(1, .95)))), 
                 col = "black", .width = 0.95, point_size = 3, slab_color = "cornflowerblue") +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
    theme_bw() + 
    theme(legend.position = "none", plot.background = element_rect(fill = "transparent", color = NA)))

# Plot data of temperature-increased L in water-breathers:

(plot_wact <- ggplot(w_act, aes(x = log10(L), y = b, group = experiment2)) +  # Set up canvas with outcome variable on y-axis
    geom_line(color="gray60")+
    geom_point(aes(color = temp), size=3, stroke = 1.5, alpha = 0.5) + # Plot the actual points
    scale_color_gradient2(low = "blue", mid="lightblue", high = "red", midpoint=17, 
                          limits = c(min(all_temp$temp), max(all_temp$temp))) + # Temperature code: colour mapped here
    geom_hline(yintercept = c(2/3,1),linetype = "dashed") +
    scale_x_continuous(bquote(''~log[10]~'Metabolic level' ~italic(L)~ '('*'mg'~O[2]~ g^-1~h^-1*')'), 
                       breaks = seq(-2.5, 0.5, 1), limits = c(-2.6,0.8))+
    scale_y_continuous("Slope"  ~italic(b)~ "", breaks = seq(0.2,1.6,0.2), limits = c(0.20,1.6))+
    labs(col = expression("T ("*degree*C*")")) +
    theme_classic() +
    theme(#legend.position = "none",
      axis.title = element_text(size = 15),
      axis.ticks.length = unit(.25, "cm"),
      axis.text = element_text(size = 13),
      axis.line = element_line(size = 0.75),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()))


# altogether:

(Fig_2C <- plot_wact + theme(text = element_text(size=13),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank()) +
    annotation_custom(ggplotGrob(wa_logL_estimate), xmin = -0.55, xmax = 0.95,  
                      ymin = 1.20, 
                      ymax = 1.67))

# Plot air-breathing species data and mean estimate for log10 temperature-increased L #

# Mean estimate of log10-L effect and posterior distribution:

(aa_logL_estimate <- 
    post_act %>%
    transmute(gamma_nowater    = b_log10_L) %>%
    gather(key, value) %>%
    
    ggplot(aes(x = value, group = key, color = key, fill = key)) +
    scale_x_continuous(expression(beta[italic(L)]), expand = c(0, 0), limits = c(-0.44, 0.447)) +
    scale_y_continuous(NULL, breaks = NULL) +
    stat_halfeye(aes(fill = stat(cut_cdf_qi(cdf,.width = c(1, .95)))), 
                 col = "black", .width = 0.95, point_size = 3, slab_color = "chartreuse3") +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_brewer(palette = "Greens", direction = -1, na.translate = FALSE) +
    theme_bw() + 
    theme(legend.position = "none",plot.background = element_rect(fill = "transparent", color = NA)))

# Plot data of temperature-increased L in water-breathers:

(plot_aact <- ggplot(a_act, aes(x = log10(L), y = b, group= experiment2)) +  # Set up canvas with outcome variable on y-axis
    geom_line(color="gray60")+
    geom_point(aes(color = temp), size = 3, stroke = 1.5, alpha = 0.5) + # Plot the actual points
    scale_color_gradient2(low = "blue", mid = "lightblue", high = "red", midpoint = 17, 
                          limits = c(min(all_temp$temp), max(all_temp$temp))) + # Temperature code: colour mapped here
    geom_hline(yintercept = c(2/3,1),linetype = "dashed") +
    scale_x_continuous(bquote(''~log[10]~'Metabolic level' ~italic(L)~ '('*'mg'~O[2]~ g^-1~h^-1*')'), 
                       breaks=seq(-2.5, 0.5, 1), limits=c(-2.6,0.8))+
    scale_y_continuous("Slope"  ~italic(b)~ "", breaks=seq(0.2,1.6,0.2), limits=c(0.20,1.6))+
    labs(col = expression("T ("*degree*C*")")) +
    theme_classic() +
    theme(#legend.position = "none",
      axis.title = element_text(size = 15),
      axis.ticks.length = unit(.25, "cm"),
      axis.text = element_text(size = 13),
      axis.line = element_line(size = 0.75),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()))

# altogether:
(Fig_2D <- plot_aact + theme(text = element_text(size = 13),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank()) +
    annotation_custom(ggplotGrob(aa_logL_estimate), xmin = -0.55, xmax = 0.95,  
                      ymin = 1.20, 
                      ymax = 1.67))

# Merge all plots into a final figure:
(Fig_2AB.p <- plot_grid(
  Fig_2A + theme(legend.position = "none", plot.title = element_blank(), axis.title.x=element_blank()),
  Fig_2B + theme(legend.position = "none", plot.title = element_blank()),
  align = 'v',
  labels = c("A", "B"),
  label_size = 18,
  hjust = -5,
  ncol = 1,
  rel_heights = c(0.49, 0.51)))

(Fig_2CD.p <- plot_grid(
  Fig_2C + theme(legend.position = "none", plot.title = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()),
  Fig_2D + theme(legend.position = "none", plot.title = element_blank(), axis.title.y=element_blank()),
  align = 'v',
  labels = c("C", "D"),
  label_size = 18,
  ncol = 1,
  hjust = -4,
  rel_heights = c(0.49, 0.51)))

# get legend (Temperature-coded colour)
legend_Fig_2 <- get_legend(
  # create some space to the left of the legend
  plot_wtemp + theme(legend.title = element_text(size = 14),
                     legend.text = element_text(size = 12),
                     legend.box.margin = margin(0, 0, 0,0.5))
)

(Fig_2_nolegend <- plot_grid(
  Fig_2AB.p,
  Fig_2CD.p,
  align = 'h',
  rel_widths = c(0.52, 0.48)))

(Fig_2 <- plot_grid(Fig_2_nolegend, legend_Fig_2, rel_widths = c(0.93,0.07)))

# ------------------------------------------------------------------------------

# Figure 3 ####

# Graphical representation of the quantitative prediction 
# on the change in slopes b raising with activity-increased L 

# Generate a function to simulate the increase in b with activity, from b at resting MR to b at MMR:

simulatedata <- function(n, massrange, act_level, L_rmr, FAS, b_rmr, b_mmr) {
  
  # Generate a log sequence of numbers for equal-spaced values 
  # of the mass range in a log scale:
  log_M <- seq(from = log10(massrange[1]), to = log10(massrange[2]), length.out = n)
  
  # Generate (anti-log) mass values:
  x <- 10^(log_M)
  mass <- rep(x, length(act_level)) # same length as activity level values
  
  # Model activity level:
  Act_L <- sort(rep(act_level, n))
  
  # Calculate the geometric mean of the mass range:
  Mm = exp((log(min(mass)) + log(max(mass)))/2)
  
  # Calculate the scaling coefficient (a) 
  # of resting (non-active) metabolic rate (rmr)
  # given a reference value of L at rest (L rmr):
  a_rmr = (L_rmr * Mm) / (Mm ^ b_rmr)
  
  # Resting metabolic rate at a reference L at rest:
  R_rmr = a_rmr * mass ^ b_rmr
  
  # Model metabolic level at maximal metabolic rate (L mmr),
  # calculated as L at rest multiplied by a factorial aerobic scope (FAS)
  # i.e., ratio between MMR / RMR:
  L_mmr = L_rmr * FAS
  
  # Calculate aerobic scope (AS) as difference between L mmr and L rmr:
  L_AS = L_mmr - L_rmr
  
  # Model active metabolic level (L act), as the proportion of aerobic scope
  # being used by the animal, as a value of activity level (Act_L):
  L_act = Act_L * L_AS
  
  # Calculate scaling coefficient (a) of active metabolic rate (R act) 
  # at a reference value of active metabolic level (L act)
  a_act = (L_act * Mm) / (Mm ^ b_mmr)
  
  # Calculate active metabolic rate:
  R_act <- a_act * mass ^ b_mmr
  
  # Model total respiration rate as the sum of resting and
  # active metabolic rates:
  resp <- R_rmr + R_act
  
  return(data.frame(mass = mass, Act_L = Act_L, resp = resp))
}


# Give real values to the function:
# (Note that all values are adjusted at 20 degrees Celsius
# to compare between water- and air-breathing species)

# Correct metabolic levels to a common temperature of T = 20 degrees C 

# Linear mixed-effects models (LMM) for the effect of temperature on metabolic level L
# using the data from inactive or minimal active animals only:

# LMM for water-breathers:
(lmer_L_wtemp <- 
    summary(lmer(log10(L) ~ temp + (temp|experiment), 
                 control = lmerControl(optimizer ="Nelder_Mead"), data = w_temp)))

# LMM for air-breathers:
(lmer_L_atemp <- 
    summary(lmer(log10(L) ~ temp + (temp|experiment),
                 control = lmerControl(optimizer ="Nelder_Mead"), data = a_temp)))

# Calculating Q10 following estimates from the LMM

# Q10 in water-breathers:
(L_at0w <- 10^((coef(lmer_L_wtemp)[1] + coef(lmer_L_wtemp)[2]*(0)))) # estimated L at T = 0 degrees C
(L_at30w <- 10^((coef(lmer_L_wtemp)[1] + coef(lmer_L_wtemp)[2]*(30))))#  estimated L at T = 30 degrees C
#
(Q10_w <- (L_at30w/L_at0w)^(10/(30-0))) # Q10 value for water-breathing species

# Q10 in air-breathers:
(L_at0a <- 10^((coef(lmer_L_atemp)[1] + coef(lmer_L_atemp)[2]*(0)))) # estimated L at T = 0 degrees C
(L_at30a <- 10^((coef(lmer_L_atemp)[1] + coef(lmer_L_atemp)[2]*(30))))#  estimated L at T = 30 degrees C
#
(Q10_a <- (L_at30a/L_at0a)^(10/(30-0))) # Q10 value for air-breathing species

# Water-breathing species:

# get data on experiments that measured resting (minimum activity) to maximum L:
w_act_minL_maxL <- w_act %>% 
  filter(comp_L == "rest_max")

# Correct metabolic level L  (in mg O2 g-1 h-1) at 20 degrees C:
(L_at20_w <- (w_act_minL_maxL$L) * Q10_w ^ ((20 - w_act_minL_maxL$temp) / 10))

(w_min_to_max_L_Tadj <- cbind(w_act_minL_maxL, L_at20_w))

# Separate dataset between maximal L and minimal L:
w_act_maxL <- w_min_to_max_L_Tadj %>% 
  filter(L_state == "maximal")
w_act_minL <- w_min_to_max_L_Tadj %>% 
  filter(L_state == "minimal")

# Mean values of experiments for each species:

# Maximum metabolic level:
(w_act_maxL_mean_exp <- ddply(w_act_maxL, .(species, experiment2), summarize,
                              b = mean(b),
                              L = mean(L_at20_w)))

# Resting (minimum) metabolic level:
(w_act_minL_mean_exp <- ddply(w_act_minL, .(species, experiment2), summarize,
                              b = mean(b),
                              L = mean(L_at20_w)))

# Mean values of species, using mean values from experiments,
# if more than 1 experiment were available for a species:

# Max. L
(w_act_maxL_mean_spp <- ddply(w_act_maxL_mean_exp, .(species), summarize,
                              b = mean(b),
                              L = mean(L)))

# Min. L:
(w_act_minL_mean_spp <- ddply(w_act_minL_mean_exp, .(species), summarize,
                              b = mean(b),
                              L = mean(L)))


# Create a dataset with mean values of max. and min. L of species:
(w_20C_spp_minL <- data.frame(state = "minimal", w_act_minL_mean_spp))
(w_20C_spp_maxL <- data.frame(state = "maximal", w_act_maxL_mean_spp))   

(w_20C_spp <- rbind(w_20C_spp_maxL,w_20C_spp_minL)) # water-breathers' dataset

# Calculating overall means and 95% confidence intervals for min. L and max. L:

# Metabolic level (L)
(w_min.L <- mean_sd(log10(w_act_minL_mean_spp$L))) # in log10
(w_max.L <- mean_sd(log10(w_act_maxL_mean_spp$L))) # in log10

## Slope b
(w_b.min <- mean_sd(w_act_minL_mean_spp$b))
(w_b.max <- mean_sd(w_act_maxL_mean_spp$b))

# Dataset of mean b and L values:
(df_bL_20C_w <- data.frame(state = c("minimal","maximal"),
                           b_mean = c(w_b.min$y, w_b.max$y),
                           b_max = c(w_b.min$ymax, w_b.max$ymax),
                           b_min = c(w_b.min$ymin, w_b.max$ymin),
                           log10_L_mean = c(w_min.L$y, w_max.L$y),
                           log10_L_max = c(w_min.L$ymax, w_max.L$ymax),
                           log10_L_min = c(w_min.L$ymin, w_max.L$ymin)))


# Quantification of the relationship between b and L (derived from Glazier 2008, 2009, 2010)
mmid <- exp((log(1000) + log(0.1)) / 2)

# Mean slope b of resting (minimal) metabolic level L among species:
b_min_w = mean(w_b.min$y)

# Mean resting (minimal) metabolic level (L_min) among species: 
L_min_w = 10^(w_min.L$y) # in mg O2 g-1 h-1

# (Theoretical) maximal metabolic level (L_max):
L_max_w = L_min_w * 10

# Calculate scaling coefficient (a) to get inactive metabolic rate (R_min)
# at the geometric mass midpoint (mmid):
a_w <- (L_min_w * mmid) / (mmid ^ b_min_w)

# Calculate maximal a' (a'_max) to get max. metabolic rate 
aprime_max_w <- L_max_w - L_min_w

# Generate a sequence of 100 values to recreate the increase from
# R_min to max. metabolic rate with activity (R_act):
aprime_w <- seq(from = 0, to = aprime_max_w, length.out = 100)

# Obtain metabolic level (L) values:
L_w <- a_w * mmid ^ (b_min_w - 1) + aprime_w

# Calculate R_min:
R_min_w <- a_w * mmid ^ b_min_w

# Calculate total metabolic rate (R_min + R_act):
Rtotal_w <- R_min_w + aprime_w * mmid

# Get the change in slope b as metabolic rate increases with activity

# (first derivative of log total resp with respect to log m 
# evaluated at geometric midpoint)

# b_hat_w <- 1 + (b_min_w - 1) * R_min_w / Rtotal_w # eq. [5] in Appendix  (A6)
b_hat_w <- 1 + (b_min_w - 1) * (L_min_w / L_w) # eq. [6] in Appendix  (A6)

# data frame:
df_w <- data.frame(b_hat_w, L_w, aprime_w)

plot(b_hat_w ~ log10(L_w), type = "l", xlab = "log L", 
     ylab = "dlog(total r) / dlog(m) at geometric mass midpoint",
     data = df_w)


# Air-breathing species:

# get data on experiments that measured resting (minimum activity) to maximum L:
a_act_minL_maxL <- a_act %>% 
  filter(comp_L == "rest_max")

# Correct metabolic level L  (in mg O2 g-1 h-1) at 20 degrees C:
(L_at20_a <- (a_act_minL_maxL$L) * Q10_a ^ ((20 - a_act_minL_maxL$temp) / 10))

(a_min_to_max_L_Tadj <- cbind(a_act_minL_maxL, L_at20_a))

# Separate dataset between maximal L and minimal L:
a_act_maxL <- a_min_to_max_L_Tadj %>% 
  filter(L_state == "maximal")
a_act_minL <- a_min_to_max_L_Tadj %>% 
  filter(L_state == "minimal")

# Mean values of experiments for each species:

# Maximum metabolic level:
(a_act_maxL_mean_exp <- ddply(a_act_maxL, .(species, experiment2), summarize,
                              b = mean(b),
                              L = mean(L_at20_a)))

# Resting (minimum) metabolic level:
(a_act_minL_mean_exp <- ddply(a_act_minL, .(species, experiment2), summarize,
                              b = mean(b),
                              L = mean(L_at20_a)))

# Mean values of species, using mean values from experiments,
# if more than 1 experiment were available for a species:

# Max. L
(a_act_maxL_mean_spp <- ddply(a_act_maxL_mean_exp, .(species), summarize,
                              b = mean(b),
                              L = mean(L)))

# Min. L:
(a_act_minL_mean_spp <- ddply(a_act_minL_mean_exp, .(species), summarize,
                              b = mean(b),
                              L = mean(L)))


# Create a dataset with mean values of max. and min. L of species:
(a_20C_spp_minL <- data.frame(state = "minimal", a_act_minL_mean_spp))
(a_20C_spp_maxL <- data.frame(state = "maximal", a_act_maxL_mean_spp))   

(a_20C_spp <- rbind(a_20C_spp_maxL, a_20C_spp_minL)) # water-breathers' dataset

# Calculating overall means and 95% confidence intervals for min. L and max. L:

# Metabolic level (L)
(a_min.L <- mean_sd(log10(a_act_minL_mean_spp$L))) # in log10
(a_max.L <- mean_sd(log10(a_act_maxL_mean_spp$L))) # in log10

## Slope b
(a_b.min <- mean_sd(a_act_minL_mean_spp$b))
(a_b.max <- mean_sd(a_act_maxL_mean_spp$b))

# Dataset of mean b and L values:
(df_bL_20C_a <- data.frame(state = c("minimal","maximal"),
                           b_mean = c(a_b.min$y, a_b.max$y),
                           b_max = c(a_b.min$ymax, a_b.max$ymax),
                           b_min = c(a_b.min$ymin, a_b.max$ymin),
                           log10_L_mean = c(a_min.L$y, a_max.L$y),
                           log10_L_max = c(a_min.L$ymax, a_max.L$ymax),
                           log10_L_min = c(a_min.L$ymin, a_max.L$ymin)))

# Quantification of the relationship between b and L (derived from Glazier 2008, 2009, 2010)
mmid <- exp((log(1000) + log(0.1)) / 2)

# Mean slope b of resting (minimal) metabolic level L among species:
b_min_a = mean(a_b.min$y)

# Mean resting (minimal) metabolic level (L_min) among species: 
L_min_a = 10^(a_min.L$y) # in mg O2 g-1 h-1

# (Theoretical) maximal metabolic level (L_max):
L_max_a = L_min_a * 10

# Calculate scaling coefficient (a) to get inactive metabolic rate (R_min)
# at the geometric mass midpoint (mmid):
a_a <- (L_min_a * mmid) / (mmid ^ b_min_a)

# Calculate maximal a' (a'_max) to get max. metabolic rate 
aprime_max_a <- L_max_a - L_min_a

# Generate a sequence of 100 values to recreate the increase from
# R_min to max. metabolic rate with activity (R_act):
aprime_a <- seq(from = 0, to = aprime_max_a, length.out = 100)

# Obtain metabolic level (L) values:
L_a <- a_a * mmid ^ (b_min_a - 1) + aprime_a

# Calculate R_min:
R_min_a <- a_a * mmid ^ b_min_a

# Calculate total metabolic rate (R_min + R_act):
Rtotal_a <- R_min_a + aprime_a * mmid

# Get the change in slope b as metabolic rate increases with activity

# (first derivative of log total resp with respect to log m 
# evaluated at geometric midpoint)

# b_hat_a <- 1 + (b_min_a - 1) * R_min_a / Rtotal_a # eq. [5] in Appendix  (A6)
b_hat_a <- 1 + (b_min_a - 1) * (L_min_a / L_a) # eq. [6] in Appendix  (A6)

# data frame:
df_a <- data.frame(b_hat_a, L_a, aprime_a)

# plot to check
plot(b_hat_a ~ log10(L_a), type = "l", xlab = "log L", 
     ylab = "dlog(total r) / dlog(m) at geometric mass midpoint",
     data = df_a)


# PLOTS:

# Water-breathing species:

(Fig_3A.p <- ggplot(df_bL_20C_w, aes(y= b_mean, x=10^(log10_L_mean)), group = state) +  #
    geom_line(data = w_20C_spp, aes(y = b, x = L, group = species), color="gray80", alpha = 0.5) + 
    geom_point(data = w_20C_spp, aes(y = b, x = L, group = species, col = state), size=2, alpha=0.5) +
    geom_errorbar(aes(ymax = b_max,
                      ymin = b_min, col = state), width = 0, size=1.5, alpha=0.5) +
    geom_errorbarh(aes(xmax = 10^(log10_L_max),
                       xmin = 10^(log10_L_min), col = state), height = 0, size=1.5, alpha=0.45) +
    geom_point(aes(fill=state, col =state), alpha=1, size=5, pch=21) + #
    scale_colour_manual(values=c("purple", "green")) +
    scale_fill_manual(values=c("purple", "green")) +
    geom_hline(yintercept = c(2/3,1),linetype = "dashed") +
    scale_y_continuous("Slope "~italic(b)~ "", breaks=seq(0.5, 1.2, 0.10), limits=c(0.45,1.25))+
    scale_x_log10(bquote(''~italic(L)~ '('*'mg'~O[2]~ g^-1~h^-1*') at 20°C'), limits=c(0.023,1.6)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size=13),
          plot.title = element_text(size=17),
          axis.ticks.length = unit(.25, "cm"),
          axis.text = element_text(size=13),
          axis.line = element_line(size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

# plot prediction of the raise in b with activity increased-L by 10-fold:
(Fig_3A <- 
    Fig_3A.p + 
    geom_line(data = df_w, aes(x = L_w, y = b_hat_w),
              lty = 4, alpha=0.45, size=1.2))

# Air-breathing species:

(Fig_3B.p <- ggplot(df_bL_20C_a, aes(y= b_mean, x=10^(log10_L_mean)), group = state) +  #
    geom_line(data = a_20C_spp, aes(y = b, x = L, group = species), color="gray80", alpha = 0.5) +
    geom_point(data = a_20C_spp, aes(y = b, x = L, group = species, col = state), size=2, alpha=0.5) +
    geom_errorbar(aes(ymax = b_max,
                      ymin = b_min, col = state), width = 0, size=1.5, alpha=0.45) +
    geom_errorbarh(aes(xmax = 10^(log10_L_max),
                       xmin = 10^(log10_L_min), col = state), height = 0, size=1.5, alpha=0.5) +
    geom_point(aes(fill=state, col =state), alpha=1, size=5, pch=21) + #
    scale_colour_manual(values=c("purple", "green")) + # Activity level
    scale_fill_manual(values=c("purple", "green")) + # Activity level
    geom_hline(yintercept = c(2/3,1),linetype = "dashed") +
    scale_y_continuous("Slope "~italic(b)~ "", breaks=seq(0.5, 1.2, 0.10), limits=c(0.45,1.25))+
    scale_x_log10(bquote(''~italic(L)~ '('*'mg'~O[2]~ g^-1~h^-1*') at 20°C'), limits=c(0.023,1.6)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size=13),
          plot.title = element_text(size=17),
          axis.ticks.length = unit(.25, "cm"),
          axis.text = element_text(size=13),
          axis.line = element_line(size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

# plot prediction of the raise in b with activity increased-L by 10-fold:
(Fig_3B <- 
    Fig_3B.p + 
    geom_line(data = df_a, aes(x = L_a, y = b_hat_a),
              lty = 4, alpha=0.45, size=1.2))

# Merge plots for water- and air-breathing species:
(Fig_3 <- plot_grid(
  Fig_3A,
  Fig_3B + theme(axis.title.y = element_blank()),
  align = 'h',
  labels = c("A", "B"),
  nrow = 1,
  rel_widths = c(0.49, 0.51)))

# ------------------------------------------------------------------------------

# Figure 4 ####

# Get mean estimates and posterior distributions of the effect of log10 L
# from brms models:

# Effect of temperature-increased L:
(temp_effect_L <- 
   post_temp %>%
   transmute(intercept = b_Intercept,
             gamma_water    = b_log10_L + `b_log10_L:groupwaterMbreather`,
             gamma_notwater = b_log10_L) %>%
   tidyr::gather(key, value) %>%
   group_by(key) %>%
   dplyr::summarise(mean = mean(value),
                    CI_lower = quantile(value, probs = c(0.05)),
                    CI_upper = quantile(value, probs = c(0.95)))%>%
   as.data.frame())

# Effect of activity-increased L:
(act_effect_L <- 
    post_act %>%
    transmute(intercept = b_Intercept,
              gamma_water    = b_log10_L + `b_log10_L:groupwaterMbreather`,
              gamma_notwater = b_log10_L) %>%
    tidyr::gather(key, value) %>%
    group_by(key) %>%
    dplyr::summarise(mean = mean(value),
                     CI_lower = quantile(value, probs = c(0.05)),
                     CI_upper = quantile(value, probs = c(0.95))) %>%
    as.data.frame())


# Effect of temperature-increased L ####

# Data preparation

# Convert all metabolic rates in mg O2 h-1:
# Using metabolic levels (metabolic rate at the the geometric mass midpoint
# of each regression, in mg O2 g-1 h-1), we can calculate the scaling coefficient
# (a), in mg O2 g-b h-1 for each regression, and then use the slopes b to estimate 
# the metabolic rates at the minimum and maximum body mass from those regressions

# Water-breathers:
(a_mg.g.h_w.temp <- (w_temp$L * w_temp$mass_mid.g) / ((w_temp$mass_mid.g) ^ w_temp$b))
(Rmin_w.temp <- a_mg.g.h_w.temp * (w_temp$min_mass.g) ^ w_temp$b)
(Rmax_w.temp <- a_mg.g.h_w.temp * (w_temp$max_mass.g) ^ w_temp$b)

# Air-breathers:
(a_mg.g.h_a.temp <- (a_temp$L * a_temp$mass_mid.g) / ((a_temp$mass_mid.g) ^ a_temp$b))
(Rmin_a.temp <- a_mg.g.h_a.temp * (a_temp$min_mass.g) ^ a_temp$b)
(Rmax_a.temp <- a_mg.g.h_a.temp * (a_temp$max_mass.g) ^ a_temp$b)

# Calculate upper and lower boundaries for the slope b at minimum and maximum 
# temperature-increased L, according to the model's mean estimates:

# Water-breathers:
(b_Lmin.wtemp <- as.numeric(temp_effect_L[3,2] + temp_effect_L[2,2] * log10(min(w_temp$L))))
(b_Lmax.wtemp <- as.numeric(temp_effect_L[3,2] + temp_effect_L[2,2] * log10(max(w_temp$L))))

# The estimated change in b with a 10-fold increase in L by warming:
(b_Lmin.wtemp_L0.1 <- as.numeric(temp_effect_L[3,2] + temp_effect_L[2,2] * log10(0.01)))
(b_Lmin.wtemp_L1 <- as.numeric(temp_effect_L[3,2] + temp_effect_L[2,2] * log10(0.1)))
# b_Lmin.wtemp_L0.1/b_Lmin.wtemp_L1
b_Lmin.wtemp_L1 - b_Lmin.wtemp_L0.1

# Air-breathers:

# We found no consistent change in b with log10 L:
(b_Lmin.atemp <- as.numeric(temp_effect_L[3,2] + temp_effect_L[1,2] * log10(min(a_temp$L))))
(b_Lmax.atemp <- as.numeric(temp_effect_L[3,2] + temp_effect_L[1,2] * log10(max(a_temp$L))))

# mean b, grouping by experiments, then by species:
spp_mean_a <- 
  a_temp %>%
  group_by(experiment) %>%
  group_by(species_phylo) %>%
  dplyr::summarise(mean_b = mean(b), mean_log10L = mean(log10_L))

mean_sd(spp_mean_a$mean_b) # overall mean b and standard deviation within species
(mean_b_a <- mean(spp_mean_a$mean_b))
(mean_log10L_a <- mean(spp_mean_a$mean_log10L))

# The estimated change in b with a 10-fold increase in L by warming:
# (b_Lmin.atemp_L0.1 <- as.numeric(temp_effect_L[3,2] + temp_effect_L[1,2] * log10(0.01)))
# (b_Lmin.atemp_L1 <- as.numeric(temp_effect_L[3,2] + temp_effect_L[1,2] * log10(0.1)))
# b_Lmin.atemp_L0.1/b_Lmin.atemp_L1

# Calculate boundaries at a reference mass:
# the midpoint of the mass range covered in each dataset (water- or air-breathing species)
(mass_mid_w.temp <- exp((log((min(w_temp$min_mass.g)))+log(max(w_temp$max_mass.g)))/2)) # water-breathing spp.
(mass_mid_a.temp <- exp((log((min(a_temp$min_mass.g)))+log(max(a_temp$max_mass.g)))/2)) # air-breathing spp.

# calculate scaling coefficient (a) for boundaries, 
# at max. and min. metabolic level L of each dataset:

# Water-breathers:
(a_minL_wtemp <- as.numeric((min(w_temp$L) * mass_mid_w.temp)/(mass_mid_w.temp ^ b_Lmin.wtemp)))
(a_maxL_wtemp <- as.numeric((max(w_temp$L) * mass_mid_w.temp)/(mass_mid_w.temp ^ b_Lmax.wtemp)))

# Air-breathers:
(a_minL_atemp <- as.numeric((min(a_temp$L) * mass_mid_a.temp)/(mass_mid_a.temp ^ b_Lmin.atemp)))
(a_maxL_atemp <- as.numeric((max(a_temp$L) * mass_mid_a.temp)/(mass_mid_a.temp ^ b_Lmax.atemp)))
(a_meanL_atemp <- as.numeric((10^(mean_log10L_a) * mass_mid_a.temp)/(mass_mid_a.temp ^ mean(mean_b_a))))

# Effect of activity-increased L ####

# Data preparation

# Correct metabolic levels to mean temperature of each dataset:
summary(w_act$temp) # 15 degrees C for water-breathers' data
summary(a_act$temp) # 25 degrees C for air-breathers' data

# Water-breathers:
(L_at15_w <- (w_act$L) * Q10_w ^ ((15 - w_act$temp) / 10)) # L (in mg O2 g-1 h-1) at 15 degrees C

# Metabolic rate (in mg O2 h-1) at 15 degrees C, estimated at the mass-midpoint of each regression:
(R_at15_w <- L_at15_w * (w_act$mass_mid.g))

# Air-breathers:
(L_at25_a <- (a_act$L) * Q10_w ^ ((15 - a_act$temp) / 10)) # L (in mg O2 g-1 h-1) at 25 degrees C

# Metabolic rate (in mg O2 h-1) at 25 degrees C, estimated at the mass-midpoint of each regression:
(R_at25_a <- L_at25_a * (a_act$mass_mid.g))

# Convert all metabolic rates in mg O2 h-1 corrected to a common temperature:
# Using metabolic levels (metabolic rate at the the geometric mass midpoint
# of each regression, in mg O2 g-1 h-1), we can calculate the scaling coefficient
# (a), in mg O2 g-b h-1 for each regression, and then use the slopes b to estimate 
# the metabolic rates at the minimum and maximum body mass from those regressions

# Water-breathers:
(a_mg.g.h_w.act <- (R_at15_w) / ((w_act$mass_mid.g) ^ w_act$b))
(Rmin_w.act <- a_mg.g.h_w.act * (w_act$min_mass.g) ^ w_act$b)
(Rmax_w.act <- a_mg.g.h_w.act * (w_act$max_mass.g) ^ w_act$b)

# Air-breathers:
(a_mg.g.h_a.act <- (R_at25_a) / ((a_act$mass_mid.g) ^ a_act$b))
(Rmin_a.act <- a_mg.g.h_a.act * (a_act$min_mass.g) ^ a_act$b)
(Rmax_a.act <- a_mg.g.h_a.act * (a_act$max_mass.g) ^ a_act$b)

# Calculate upper and lower boundaries for the slope b at minimum and maximum 
# temperature-increased L, according to the model's mean estimates:

# Water-breathers:
(b_Lmin.wact <- as.numeric(act_effect_L[3,2] + act_effect_L[2,2] * log10(min(L_at15_w))))
(b_Lmax.wact <- as.numeric(act_effect_L[3,2] + act_effect_L[2,2] * log10(max(L_at15_w))))

# mean b, grouping by experiments, then by species:
spp_mean_w <- 
  w_act %>%
  group_by(experiment) %>%
  group_by(species_phylo) %>%
  dplyr::summarise(mean_b = mean(b), mean_log10L = mean(log10_L))

mean_sd(spp_mean_w$mean_b) # overall mean b and standard deviation within species
(mean_b_w <- mean(spp_mean_w$mean_b))
(mean_log10L_w <- mean(spp_mean_w$mean_log10L))

# The estimated change in b with a 10-fold increase in L by warming:
# (b_L_act_L0.1w <- as.numeric(act_effect_L[3,2] + act_effect_L[2,2] * log10(0.1)))
# (b_L_act_L1w <- as.numeric(act_effect_L[3,2] + act_effect_L[2,2] * log10(1)))
# b_L_act_L1w/b_L_act_L0.1w

# Air-breathers:
(b_Lmin.aact <- as.numeric(act_effect_L[3,2] + act_effect_L[1,2] * log10(min(L_at25_a))))
(b_Lmax.aact <- as.numeric(act_effect_L[3,2] + act_effect_L[1,2] * log10(max(L_at25_a))))

# The estimated change in b with a 10-fold increase in L by warming:
(b_L_act_L0.01a <- as.numeric(act_effect_L[3,2] + act_effect_L[1,2] * log10(0.01)))
(b_L_act_L0.1a <- as.numeric(act_effect_L[3,2] + act_effect_L[1,2] * log10(0.1)))
# b_L_act_L1a/b_L_act_L0.1a
b_L_act_L0.1a - b_L_act_L0.01a

# Calculate boundaries at a reference mass:
# the midpoint of the mass range covered in each dataset (water- or air-breathing species)
(mass_mid_w.act <- exp((log(min(w_act$min_mass.g)) + log(max(w_act$max_mass.g))) / 2))
(mass_mid_a.act <- exp((log(min(a_act$min_mass.g)) + log(max(a_act$max_mass.g))) / 2))

# calculate scaling coefficient (a) for boundaries, 
# at max. and min. metabolic level L of each dataset:

# Water-breathers:
(a_minL_wact <- as.numeric((min(L_at15_w) * mass_mid_w.act)/(mass_mid_w.act ^ b_Lmin.wact)))
(a_maxL_wact <- as.numeric((max(L_at15_w) * mass_mid_w.act)/(mass_mid_w.act ^ b_Lmax.wact)))
(a_meanL_wact <- as.numeric((10^(mean_log10L_w) * mass_mid_w.act)/(mass_mid_w.act ^ mean_b_w)))

# Air-breathers:
(a_minL_aact <- as.numeric((min(L_at25_a) * mass_mid_a.act)/(mass_mid_a.act ^ b_Lmin.aact)))
(a_maxL_aact <- as.numeric((max(L_at25_a) * mass_mid_a.act)/(mass_mid_a.act ^ b_Lmax.aact)))

# PLOT: effect temperature-increased metabolic level (L) on b:

# Water-breathing species data:
(Fig_4A.p <-
    ggplot(w_temp,aes(y=log10(L),x=log10(mass_mid.g))) +
    geom_segment(aes(x = log10(min_mass.g), 
                     xend = log10(max_mass.g), 
                     y = log10(Rmin_w.temp), 
                     yend = log10(Rmax_w.temp),
                     col = temp), 
                 size = 0.71,
                 alpha = 0.6) +
    scale_color_gradient2(low = "blue", mid="lightblue", high = "red", midpoint=17, 
                          limits = c(min(all_temp$temp), max(all_temp$temp))) + # Temperature code: colour mapped here
    geom_abline(intercept = log10(a_minL_wtemp), slope = b_Lmin.wtemp, col="black", lty=2, lwd=1.1, alpha=0.62)+
    geom_abline(intercept = log10(a_maxL_wtemp), slope = b_Lmax.wtemp, col="black", lty=2, lwd=1.1, alpha=0.62)+ 
    scale_x_continuous(bquote(''~log[10]~'Body mass (g)'), 
                       breaks=seq(-2, 4, 2), limits=c(-2.8,5.2))+
    scale_y_continuous(bquote(''~log[10]~'Metabolic rate ('*' mg'~O[2]~h^-1*')'), 
                       breaks=seq(-6, 4, 2), limits=c(-6,4.3))+
    theme_classic() +
    theme(axis.title = element_text(size=14),
          axis.ticks.length = unit(.25, "cm"),
          axis.text = element_text(size=13),
          axis.line = element_line(size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

# Air-breathing species data:
(Fig_4B.p <-
    ggplot(a_temp,aes(y=log10(L),x=log10(mass_mid.g))) +
    geom_segment(aes(x = log10(min_mass.g), 
                     xend = log10(max_mass.g), 
                     y = log10(Rmin_a.temp), 
                     yend = log10(Rmax_a.temp),
                     col = temp), 
                 size = 0.71,
                 alpha = 0.6) +
    scale_color_gradient2(low = "blue", mid="lightblue", high = "red", midpoint=17, 
                          limits = c(min(all_temp$temp), max(all_temp$temp))) + # Temperature code: colour mapped here
    #geom_abline(intercept = log10(a_minL_atemp), slope = b_Lmin.atemp, col="black", lty=2, lwd=1.1, alpha=0.62)+ 
    #geom_abline(intercept = log10(a_maxL_atemp), slope = b_Lmax.atemp, col="black", lty=2, lwd=1.1, alpha=0.62)+
    geom_abline(intercept = log10(a_meanL_atemp), slope = mean_b_a, col="black", lty=2, lwd=1.1, alpha=0.62)+
    scale_x_continuous(bquote(''~log[10]~'Body mass (g)'), 
                       breaks=seq(-2, 4, 2), limits=c(-2.8,5.2))+
    scale_y_continuous(bquote(''~log[10]~'Metabolic rate ('*' mg'~O[2]~h^-1*')'), 
                       breaks=seq(-6, 4, 2), limits=c(-6,4.3))+
    labs(col = expression("T ("*degree*C*")")) +
    theme_classic() +
    theme(axis.title = element_text(size=14),
          axis.ticks.length = unit(.25, "cm"),
          axis.text = element_text(size=13),
          axis.line = element_line(size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 11),
          legend.box.margin = margin(0, 0, 0,0.5)))


# PLOT: effect activity-increased metabolic level (L) on b:

# Water-breathing species data:
(Fig_4C.p <-
    ggplot(w_act,aes(y=log10(R_at15_w),x=log10(mass_mid.g))) +
    geom_segment(aes(x = log10(min_mass.g), 
                     xend = log10(max_mass.g), 
                     y = log10(Rmin_w.act), 
                     yend = log10(Rmax_w.act)), 
                 size = 0.71,
                 alpha = 0.5,
                 col= "grey70") +
    #geom_abline(intercept = log10(a_minL_wact), slope = b_Lmin.wact, col="black", lty=2, lwd=1.1, alpha=0.65)+ #
    #geom_abline(intercept = log10(a_maxL_wact), slope = b_Lmax.wact, col="black", lty=2, lwd=1.1, alpha=0.65)+ #
    # mean activity-increased log10 L:
    geom_abline(intercept = log10(a_meanL_wact), slope = mean_b_w, col="black", lty=2, lwd=1.1, alpha=0.65)+ # 
    scale_x_continuous(bquote(''~log[10]~'Body mass (g)'), 
                       breaks=seq(-4, 4, 2), limits=c(-4.5,5.2))+
    scale_y_continuous(bquote(''~log[10]~'Metabolic rate ('*' mg'~O[2]~h^-1*') at 15°C'), 
                       breaks=seq(-6, 4, 2), limits=c(-6,4.3))+
    theme_classic() +
    theme(axis.title = element_text(size=14),
          axis.ticks.length = unit(.25, "cm"),
          axis.text = element_text(size=13),
          axis.line = element_line(size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))


# Air-breathing species data:
(Fig_4D.p <-
    ggplot(a_act,aes(y=log10(R_at25_a),x=log10(mass_mid.g))) +
    geom_segment(aes(x = log10(min_mass.g), 
                     xend = log10(max_mass.g), 
                     y = log10(Rmin_a.act), 
                     yend = log10(Rmax_a.act)), 
                 size = 0.75,
                 alpha = 0.7,
                 col= "grey70") +
    geom_abline(intercept = log10(a_minL_aact), slope = b_Lmin.aact, col="black", lty=2, lwd=1.1, alpha=0.62)+ # experiment
    geom_abline(intercept = log10(a_maxL_aact), slope = b_Lmax.aact, col="black", lty=2, lwd=1.1, alpha=0.62)+ # experiment
    scale_x_continuous(bquote(''~log[10]~'Body mass (g)'), 
                       breaks=seq(-4, 4, 2), limits=c(-4.5,5.2))+
    scale_y_continuous(bquote(''~log[10]~'Metabolic rate ('*' mg'~O[2]~h^-1*') at 25°C'), 
                       breaks=seq(-6, 4, 2), limits=c(-6,4.3))+
    theme_classic() +
    theme(axis.title = element_text(size=14),
          axis.ticks.length = unit(.25, "cm"),
          axis.text = element_text(size=13),
          axis.line = element_line(size = 0.75),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))


# Merge plots into final figure

# Plots for temperature effect (A and B):
(Fig_4AB.p <- plot_grid(
  Fig_4A.p + theme(axis.title.x = element_blank(), plot.title = element_blank(), legend.position = "none"),
  Fig_4B.p + theme(plot.title = element_blank(), legend.position = c(0.9,0.25)),
  align = 'v',
  labels = c("A", "B"),
  label_size = 18,
  hjust = -1,
  ncol = 1,
  rel_heights = c(0.5, 0.5)))

# add title:
(title_temp <- ggdraw() + 
    draw_label(
      "Temperature",
      fontface = 'bold',
      size = 17,
      x = 0.2) +
    theme(plot.margin = margin(0, 0, 0, 0)))

(Fig_4AB.pt <- plot_grid(
  title_temp,
  Fig_4AB.p,
  align = 'v',
  ncol = 1,
  rel_heights = c(0.05, 0.95)))

# Plots for activity effect (C and D):
(Fig_4CD.p <- plot_grid(
  Fig_4C.p + theme(plot.title = element_blank(), axis.title.x = element_blank()),
  Fig_4D.p + theme(plot.title = element_blank()),
  align = 'v',
  labels = c("C", "D"),
  label_size = 18,
  hjust = -1,
  ncol = 1,
  rel_widths = c(0.5, 0.5)))

# add title:
(title_act <- ggdraw() + 
    draw_label(
      "Activity",
      fontface = 'bold',
      size = 17,
      x = 0.165) +
    theme(plot.margin = margin(0, 0, 0, 0)))

(Fig_4CD.pt <- plot_grid(
  title_act,
  Fig_4CD.p,
  align = 'v',
  ncol = 1,
  rel_heights = c(0.05, 0.95)))


# Final figure 4:
(Fig_4 <- plot_grid(
  Fig_4AB.pt,
  Fig_4CD.pt,
  align = 'h',
  hjust = 1,
  nrow = 1,
  rel_widths = c(0.5, 0.5)))

# ------------------------------------------------------------------------------
# Save figures #####
# In pdf format (png optional) #
dev.off()

# Figure 1:

pdf(file = "../Figures/Figure 1.pdf", width = 10, height = 7.5, useDingbats = FALSE)
# png(file = "../Figures/Figure 1.png", width = 10, height = 7.5, units = "in", res = 600)

Fig_1 

dev.off()

# Figure 2:

pdf(file = "../Figures/Figure 2.pdf", width = 16, height = 11, useDingbats = FALSE)
# png(file = "../Figures/Figure 2.png",width = 16, height = 11, units = "in", res = 600)

Fig_2

dev.off()

# Figure 3:

pdf(file = "../Figures/Figure 3.pdf", width = 10, height = 5, useDingbats = FALSE)
# png(file = "../Figures/Figure 3.png", width = 10, height = 5, units = "in", res = 600)

Fig_3

dev.off()

# Figure 4:

pdf(file = "../Figures/Figure 4.pdf", width = 13, height = 10, useDingbats = FALSE)
# png(file = "../Figures/Figure 4.png", width = 13, height = 10, units = "in", res = 600)

Fig_4

dev.off()

# ------------------------------------------------------------------------------
# Save data on R session and packages versions for reproducibility shake:
sink("../R_session/Figures_ms_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################