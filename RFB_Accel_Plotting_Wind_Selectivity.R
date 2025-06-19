rm(list=ls())

#### Create plot of the outbound and inbound wind selectivity results ####

# Ruth Dunn 2025
# r.e.dunn@lancaster.ac.uk

library(tidyverse)
library(brms)      # for modelling
library(ROCR)      # for area under a curve calculations
library(pROC)      # for plotting the ROC curve
library(sf)        # for plotting Chagos shapefile
library(ncdf4)     # for loading wind data
library(patchwork) # for multi plots
library(scales)    # for adjusting axes breaks
library(raster)    # for plotting wind
library(viridis)   # for nice colour palettes

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Load data ####

dat_outbound <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Wind_Calcs.csv") %>%
  # Select outbound and at-sea points:
  filter(Trip_section == "Outbound" & atCP == "No") %>%
  # Scale and center important variables:
  mutate(sc_ws = scale(wind_speed),
         sc_wd = scale(relative_wind_dir)) %>%
  mutate(point_id = case_when(replicate == "Original" ~ "1",
                              TRUE ~ str_extract(as.character(replicate), "\\d+"))) %>%
  # Revalue 'Presence' column so that present = 1 and absent = 0
  mutate(Replicate = str_replace(as.character(replicate), "Replicate\\s?\\d+", "Replicate")) %>%
  mutate(Presence = as.numeric(dplyr::recode(Replicate, Replicate = "0", Original = "1"))) %>%
  dplyr::select(Individual, ID, Sex, DateTime, Replicate, Presence, point_id,
                Lon, Lat, wind_speed, relative_wind_dir, sc_ws, sc_wd) %>%
  na.omit()

summary(dat_outbound$relative_wind_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT INTERACTION (OUTBOUND) ####

full_mod_out <- readRDS(paste("RFB_Accel_Models/Wind_selection_i/Wind_selection_interaction_gam_4k_", 10, ".rds", sep = ""))

new_plot_dat <- expand.grid(
  sc_ws = seq(min(dat_outbound$sc_ws, na.rm = TRUE), max(dat_outbound$sc_ws, na.rm = TRUE), length.out = 6),
  sc_wd = seq(min(dat_outbound$sc_wd, na.rm = TRUE), max(dat_outbound$sc_wd, na.rm = TRUE), length.out = 20)) %>%
  mutate(
    # Add predicted mean speed
    predicted_speed = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, mean),
    # Un-scale and un-center variables
    wind_speed = (sc_ws * sd(dat_outbound$wind_speed, na.rm = TRUE)) + mean(dat_outbound$wind_speed, na.rm = TRUE),
    relative_wind_dir = (sc_wd * sd(dat_outbound$relative_wind_dir, na.rm = TRUE)) + mean(dat_outbound$relative_wind_dir, na.rm = TRUE)) %>%
  mutate(
    # Add uncertainty intervals
    lower = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975)) %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 2),  # Defines bin edges (0-2, 2-4, ..., 10-12)
                              labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12"),  # Labels for bins
                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

dat_outbound <- dat_outbound %>%
  mutate(Presence_points = ifelse(Presence == 1, 1.01,
                                  ifelse(Presence == 0, -0.01, NA))) %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 2),  # Defines bin edges (0-2, 2-4, ..., 10-12)
                              labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12"),  # Labels for bins
                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

mod_plot <- ggplot(dat_outbound) +
  geom_point(aes(y = Presence_points, x = relative_wind_dir, col = wind_speed_bin), alpha = 0.2, size = 0.75) +
  geom_ribbon(data = new_plot_dat,
              aes(x = relative_wind_dir, ymin = lower, ymax = upper, group = wind_speed_bin, fill = wind_speed_bin),
              alpha = 0.3) +
  geom_line(data = new_plot_dat, 
            aes(x = relative_wind_dir, y = predicted_speed, group = wind_speed_bin, col = wind_speed_bin), linewidth = 1) +
  scale_colour_viridis_d(option = "C", name = expression("Wind speed (ms"^-1*")"),
                         guide = guide_legend(nrow = 1)) +
  scale_fill_viridis_d(option = "C", name = expression("Wind speed (ms"^-1*")"),
                       guide = guide_legend(nrow = 1)) +
  theme_bw() +
  labs(y = expression("Probability of selection"),
       x = expression("Relative wind direction (tailwind \u2194 headwind)"),
       subtitle = "a Outbound") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        legend.margin = margin(),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0,45,90,135,180), labels = c(0,45,90,135,180), limits = c(0,180))

# mod_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT INTERACTION (INBOUND) ####

return_mod <- readRDS("RFB_Accel_Models/Wind_selection_i/Wind_selection_inbound_interaction_gam_4k_10.rds")

dat_inbound <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Inbound_Wind_Calcs.csv") %>%
  # Select outbound and at-sea points:
  filter(Trip_section == "Return" & atCP == "No") %>%
  # Scale and center important variables:
  mutate(sc_ws = scale(wind_speed),
         sc_wd = scale(relative_wind_dir)) %>%
  mutate(point_id = case_when(replicate == "Original" ~ "1",
                              TRUE ~ str_extract(as.character(replicate), "\\d+"))) %>%
  # Revalue 'Presence' column so that present = 1 and absent = 0
  mutate(Replicate = str_replace(as.character(replicate), "Replicate\\s?\\d+", "Replicate")) %>%
  mutate(Presence = as.numeric(dplyr::recode(Replicate, Replicate = "0", Original = "1"))) %>%
  dplyr::select(Individual, ID, Sex, DateTime, Replicate, Presence, point_id,
                Lon, Lat, wind_speed, relative_wind_dir, sc_ws, sc_wd) %>%
  na.omit() %>%
  na.omit() %>%
  mutate(Presence_points = ifelse(Presence == 1, 1.01,
                                  ifelse(Presence == 0, -0.01, NA))) %>%
  mutate(Locations = ifelse(Presence == 1, "Original",
                            ifelse(Presence == 0, "Simulated", NA))) %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 2),  # Defines bin edges (0-2, 2-4, ..., 10-12)
                              labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12"),  # Labels for bins
                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

new_plot_dat_in <- expand.grid(
  sc_ws = seq(min(dat_inbound$sc_ws, na.rm = TRUE), max(dat_inbound$sc_ws, na.rm = TRUE), length.out = 6),
  sc_wd = seq(min(dat_inbound$sc_wd, na.rm = TRUE), max(dat_inbound$sc_wd, na.rm = TRUE), length.out = 20)) %>%
  mutate(
    # Add predicted mean speed
    predicted_speed = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, mean),
    # Un-scale and un-center variables
    wind_speed = (sc_ws * sd(dat_inbound$wind_speed, na.rm = TRUE)) + mean(dat_inbound$wind_speed, na.rm = TRUE),
    relative_wind_dir = (sc_wd * sd(dat_inbound$relative_wind_dir, na.rm = TRUE)) + mean(dat_inbound$relative_wind_dir, na.rm = TRUE)) %>%
  mutate(
    # Add uncertainty intervals
    lower = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975)) %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 2),  # Defines bin edges (0-2, 2-4, ..., 10-12)
                              labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12"),  # Labels for bins
                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

dat_inbound <- dat_inbound %>%
  mutate(Presence_points = ifelse(Presence == 1, 1.01,
                                  ifelse(Presence == 0, -0.01, NA)))

return_mod_plot <- ggplot(dat_inbound) +
  geom_point(aes(y = Presence_points, x = relative_wind_dir, col = wind_speed_bin), alpha = 0.2, size = 0.75) +
  geom_ribbon(data = new_plot_dat_in,
              aes(x = relative_wind_dir, ymin = lower, ymax = upper, group = wind_speed_bin, fill = wind_speed_bin),
              alpha = 0.3) +
  geom_line(data = new_plot_dat_in, 
            aes(x = relative_wind_dir, y = predicted_speed, group = wind_speed_bin, col = wind_speed_bin), linewidth = 1) +
  scale_colour_viridis_d(option = "C", name = expression("Wind speed (ms"^-1*")"),
                         guide = guide_legend(nrow = 1), drop = F) +
  scale_fill_viridis_d(option = "C", name = expression("Wind speed (ms"^-1*")"),
                       guide = guide_legend(nrow = 1), drop = F) +
  theme_bw() +
  labs(y = expression("Probability of selection"),
       x = expression("Relative wind direction (tailwind \u2194 headwind)"),
       subtitle = "b Inbound") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        legend.margin = margin(),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0,45,90,135,180), labels = c(0,45,90,135,180), limits = c(0,180))

return_mod_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Layout and save plots together ####

mod_plots <- mod_plot / return_mod_plot +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom") &
  theme(plot.tag = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        title = element_text(size = 8))

ggsave("RFB_Accel_Plots/Wind_selectivity_densities_new.png", height = 6, width = 5, dpi = 600)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
