#### Model the influence of wind on flight dynamics ####

# Ruth Dunn 2025
# r.e.dunn@lancaster.ac.uk


#Set up ####

library(tidyverse)
library(data.table)
library(brms)
library(patchwork)
library(posterior)
library(ggridges)

rm(list=ls())

# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(cmdstanr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data ####

gps_dat <- left_join(fread("RFB_Accel_Data/Accel_6_flight_variables/DG_GPS_segment_flight_info.csv") %>%
                       filter(Trip_section == "Outbound" | Trip_section == "Return") %>% # Focus on outward and return only
                       filter(Flap_duration != 0 & Flap_duration != 1) %>% # Get rid of the final point sections
                       mutate(Trip = str_extract(ID, "trip\\d+")),
                     read_csv("RFB_Accel_Data/Accel_metadata.csv") %>%
                       mutate(Individual = ID) %>%
                       dplyr::select(Individual, Mass_deployment), by = "Individual") %>%
  mutate(Flight_proportion = Flap_proportion + Glide_proportion) %>%
  filter(Flight_proportion > 0.95) %>%
  # Scale and center variables:
    mutate(wind_speed_sc = scale(wind_speed, center = T, scale = T),
         relwinddir_sc = scale(relative_wind_dir, center = T, scale = T),
         bodymass_sc = scale(Mass_deployment, center = T, scale = T)) %>%
  na.omit() %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 2),  # Defines bin edges (0-2, 2-4, ..., 10-12)
                              labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12"),  # Labels for bins
                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

# Collinearity check:

# ggplot(gps_dat) +
  # geom_point(aes(wind_speed_sc, relwinddir_sc)) # ok

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modelling ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Model flight speed ####
# 
# head(gps_dat)
# 
# speed_mod <- brm(bf(speed ~ t2(wind_speed_sc, relwinddir_sc, k = 4) +
#                       (1|Individual) + (1|Trip:Individual)),
#                  control = list(adapt_delta = 0.85),
#                  chains = 4, iter = 2000, warmup = 1000,
#                  family = gaussian(), data = gps_dat)
# 
# # Save for later:
# saveRDS(speed_mod, file="RFB_Accel_Models/Flight_speed_mod_t2.rds")
# 
# # Check convergence and fit:
# pp_check(speed_mod) # posterior predictive plot
# plot(speed_mod) # trace plots
# summary(speed_mod) # ESS & R-hat
# 
# conditional_effects(speed_mod)
# 
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # Model flap proportion ####
# 
# flap_mod <- brm(bf(Flap_proportion ~ t2(wind_speed_sc, relwinddir_sc, k = 4) +
#                      (1|Individual) + (1|Trip:Individual)),
#                 control = list(adapt_delta = 0.95),
#                 chains = 4, iter = 2000, warmup = 1000,
#                 family = Beta(), data = gps_dat)
# 
# # Save for later:
# saveRDS(flap_mod, file="RFB_Accel_Models/Flap_proportion_mod_t2.rds")
# 
# # Check convergence and fit:
# pp_check(flap_mod) # posterior predictive plot
# plot(flap_mod) # trace plots
# summary(flap_mod) # ESS & R-hat
# 
# conditional_effects(flap_mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load models:

flap_mod <- readRDS("RFB_Accel_Models/Flap_proportion_mod_t2.rds")
speed_mod <- readRDS("RFB_Accel_Models/Flight_speed_mod_t2.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create new data frames, calculating predictions:

new_speed_wind <- expand.grid(
  wind_speed = seq(0, 12, length.out = 12),
  relative_wind_dir = seq(0, 180, length.out = 14)) %>%
  # Manually scale using the same parameters as in your model
  mutate(
    wind_speed_sc = (wind_speed - mean(gps_dat$wind_speed, na.rm = TRUE)) /
      sd(gps_dat$wind_speed, na.rm = TRUE),
    relwinddir_sc = (relative_wind_dir - mean(gps_dat$relative_wind_dir, na.rm = TRUE)) /
      sd(gps_dat$relative_wind_dir, na.rm = TRUE),
  bodymass_sc = mean(gps_dat$bodymass_sc),
  Trip_section = "Outbound") %>%
  mutate(
    # Add predicted mean speed
    predicted_speed = posterior_epred(speed_mod, newdata = ., re_formula = NA) %>% apply(2, mean),
    # Un-scale and un-center variables
    wind_speed = (wind_speed_sc * sd(gps_dat$wind_speed, na.rm = TRUE)) + mean(gps_dat$wind_speed, na.rm = TRUE),
    relative_wind_dir = (relwinddir_sc * sd(gps_dat$relative_wind_dir, na.rm = TRUE)) + mean(gps_dat$relative_wind_dir, na.rm = TRUE)) %>%
  mutate(
    # Add uncertainty intervals
    lower = posterior_epred(speed_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(speed_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975)) %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 1),  # Defines bin edges (0-2, 2-4, ..., 10-12)
                              labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6",
                                         "6-7", "7-8", "8-9", "9-10", "10-11", "11-12"),  # Labels for bins
                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

new_flap_wind <- expand.grid(
  wind_speed = seq(0, 12, length.out = 12),
  relative_wind_dir = seq(0, 180, length.out = 14)) %>%
  # Manually scale using the same parameters as in your model
  mutate(
    wind_speed_sc = (wind_speed - mean(gps_dat$wind_speed, na.rm = TRUE)) /
      sd(gps_dat$wind_speed, na.rm = TRUE),
    relwinddir_sc = (relative_wind_dir - mean(gps_dat$relative_wind_dir, na.rm = TRUE)) /
      sd(gps_dat$relative_wind_dir, na.rm = TRUE),
  bodymass_sc = mean(gps_dat$bodymass_sc),
  Trip_section = "Outbound") %>%
  mutate(
    # Add predicted mean speed
    predicted_speed = posterior_epred(flap_mod, newdata = ., re_formula = NA) %>% apply(2, mean),
    # Un-scale and un-center variables
    wind_speed = (wind_speed_sc * sd(gps_dat$wind_speed, na.rm = TRUE)) + mean(gps_dat$wind_speed, na.rm = TRUE),
    relative_wind_dir = (relwinddir_sc * sd(gps_dat$relative_wind_dir, na.rm = TRUE)) + mean(gps_dat$relative_wind_dir, na.rm = TRUE)) %>%
  mutate(
    # Add uncertainty intervals
    lower = posterior_epred(flap_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(flap_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975)) %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 1),  # Defines bin edges (0-2, 2-4, ..., 10-12)
                              labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6",
                                         "6-7", "7-8", "8-9", "9-10", "10-11", "11-12"),  # Labels for bins
                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make plots of variables + wind ####

speed_wind_heatmap <- ggplot(new_speed_wind,
                             aes(x = relative_wind_dir,
                                 y = wind_speed_bin,
                                 fill = predicted_speed)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B", name = expression("Ground speed (ms"^-1*")")) +
  scale_y_discrete(name = expression("Wind speed (ms"^-1*")")) +
  scale_x_continuous(name = expression("Relative wind direction (tailwind \u2194 headwind)"),
                     breaks = c(0, 45, 90, 135, 180)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 8))
speed_wind_heatmap  

flap_wind_heatmap <- ggplot(new_flap_wind,
                            aes(x = relative_wind_dir,
                                y = wind_speed_bin,
                                fill = predicted_speed)) +
  geom_tile() +
  scale_fill_viridis_c(option = "A", name = "Proportion\nflapping") +
  scale_y_discrete(name = expression("Wind speed (ms"^-1*")")) +
  scale_x_continuous(name = expression("Relative wind direction (tailwind \u2194 headwind)"),
                     breaks = c(0, 45, 90, 135, 180)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 8))

flap_wind_heatmap

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add in viloin plots ####

gps_dat <- gps_dat %>%
  mutate(x = 1)

flap_violin <- ggplot(gps_dat, aes(x = Flap_proportion)) +
  # Fill each histogram bar using the x axis category that ggplot creates
  geom_histogram(
    aes(fill = after_stat((x))), 
    binwidth = .05, boundary = 0, color = "white") +
  # Fill with the same palette as the map
  scale_fill_viridis_c(option = "A", name = "Proportion\nflapping",
                       limits = c(min(new_flap_wind$predicted_speed),
                                  max(new_flap_wind$predicted_speed)),   # Match the limits of the heatmap
                       oob = scales::squish) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  labs(x = expression("Proportion flapping"))
flap_violin

speed_violin <- ggplot(gps_dat, aes(x = speed)) + 
  geom_histogram(aes(fill = after_stat(x)),
                 binwidth = 1, boundary = 0, colour = "white") +
  scale_fill_viridis_c(option = "B", name = expression("Ground speed (ms"^-1*")"),
                       limits = c(min(new_speed_wind$predicted_speed),
                                  max(new_speed_wind$predicted_speed)),
                       oob = scales::squish) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  labs(x = expression("Ground speed (ms"^-1*")"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save plots ####


design <- "AABB
           AABB
           CCDD"

flap_wind_heatmap + speed_wind_heatmap + flap_violin + speed_violin +
  plot_layout(design = design, axes = "collect_y") &
                    theme(legend.position = "none",
                          plot.tag = element_text(size = 9),
                          axis.title = element_text(size = 9),
                          axis.text = element_text(size = 9),
                          legend.title = element_text(size = 9),
                          legend.text = element_text(size = 9)) &
  plot_annotation(tag_levels = 'a')

ggsave("RFB_Accel_Plots/Transit_wind.png", height = 5, width = 8.5, dpi = 600)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add booby image to plot ####

library(magick)

# Call back the plot
plot <- image_read("RFB_Accel_Plots/Transit_wind.png")

# Import seabird images
RFB_img <- image_read("RFB_Accel_Plots/booby_flying.png") %>%
  image_scale("700x700")

# Stack them on top of each other
plot1 <- image_composite(plot, RFB_img, offset = "+4300+100")

# And save
image_write(plot1, "RFB_Accel_Plots/Final/Fig3_flight_behaviour_RFBimage.png")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Speed summary values for text ####

# Speed ####

mean(gps_dat$speed, na.rm = T)
sd(gps_dat$speed, na.rm = T)

# Compute conditional effects
conditional_effects_speed <- conditional_effects(speed_mod, effects = "wind_speed_sc:relwinddir_sc")
mean(conditional_effects_speed$`wind_speed_sc:relwinddir_sc`$estimate__)
mean(conditional_effects_speed$`wind_speed_sc:relwinddir_sc`$lower__)
mean(conditional_effects_speed$`wind_speed_sc:relwinddir_sc`$upper__)

# Example specific values for relwinddir_sc
specific_values_speed <- c(min(conditional_effects_speed[["wind_speed_sc:relwinddir_sc"]]$relwinddir_sc),
                     median(conditional_effects_speed[["wind_speed_sc:relwinddir_sc"]]$relwinddir_sc),
                     max(conditional_effects_speed[["wind_speed_sc:relwinddir_sc"]]$relwinddir_sc))

# Extract the relevant data
conditional_effects_speed[["wind_speed_sc:relwinddir_sc"]] %>%
  filter(relwinddir_sc %in% specific_values_speed) %>%
  group_by(relwinddir_sc) %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))

# Proportion flapping summary ####

mean(gps_dat %>%
       pull(Flap_proportion), na.rm = TRUE)
sd(gps_dat %>%
     pull(Flap_proportion), na.rm = TRUE)

flap_ce <- conditional_effects(flap_mod)
mean(flap_ce$`wind_speed_sc:relwinddir_sc`$estimate__)
mean(flap_ce$`wind_speed_sc:relwinddir_sc`$lower__)
mean(flap_ce$`wind_speed_sc:relwinddir_sc`$upper__)

# Example specific values for relwinddir_sc
specific_values_flap <- c(min(flap_ce[["wind_speed_sc:relwinddir_sc"]]$relwinddir_sc),
                           median(flap_ce[["wind_speed_sc:relwinddir_sc"]]$relwinddir_sc),
                           max(flap_ce[["wind_speed_sc:relwinddir_sc"]]$relwinddir_sc))

# Extract the relevant data
flap_ce[["wind_speed_sc:relwinddir_sc"]] %>%
  filter(relwinddir_sc %in% specific_values_flap) %>%
  group_by(relwinddir_sc) %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))
