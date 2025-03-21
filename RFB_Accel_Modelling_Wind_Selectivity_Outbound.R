rm(list=ls())

#### Habitat selection GLM of presence vs rotated data for the outbound commutes ####

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

ggplot(dat_outbound %>%
         mutate(year = lubridate::year(DateTime)) %>%
         filter(year == 2019)) +
  geom_histogram(aes(relative_wind_dir), bins = 15) +
  theme_bw() +
  labs(title = 2019) +
  ggplot(dat_outbound %>%
           mutate(year = lubridate::year(DateTime)) %>%
           filter(year == 2022)) +
  geom_histogram(aes(relative_wind_dir), bins = 15) +
  theme_bw() +
  labs(title = 2022)

# ggsave("RFB_Accel_Plots/Rel_wind_dir_hists.png")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Run full model ####

# Set up parallel processing:
# install.packages("parallelly")
options(mc.cores = parallelly::availableCores())

auc.vals <- NULL

for(i in c(1:20)){
  
  # i = 10
  
  # Subset data to include a certain number of points and weight the points by this value
  dati <- dat_outbound %>%
    filter(Presence == 1 | Presence == 0 & point_id <= i) %>%
    # Add weighting for presence/psuedo-absence:
    mutate(wt = if_else(Presence == 1, i, 1))
  
  # Fit the model:
  
  full_modi <- brm(Presence | weights(wt) ~ 0 + t2(sc_ws, sc_wd, k = 4) +
                  (1|Individual) + (1|ID:Individual),
                data = dati, family = bernoulli,
                chains = 4,
                threads = threading(4),
                backend = "cmdstanr") # use CmdStan (for speed)
  
  # Model evaluation - 'area under curve' (AUC)
  # Tutorial: https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
  AUC <- performance(prediction(predict(full_modi, type = "response")[,1],
                                as.vector(pull(na.omit(dati), Presence))),
                     measure = "auc")
  AUC <- AUC@y.values[[1]]
  # A value of 0.50 means that the model does not classify better than chance
  # A good model should have an AUC score much higher than 0.50 (preferably higher than 0.80)
  
  # Also extract the R^2 value and the number of points used:
  R2 <- bayes_R2(full_modi)[,1]
  
  # Save it up:
  auc.val <- cbind(i, AUC, R2)
  auc.vals <- rbind(auc.vals, auc.val)
  
  write_csv(as.data.frame(auc.vals), ("RFB_Accel_Models/Wind_selectivity_random_point_AUC_4k.csv"))

  # Save model:
  saveRDS(full_modi, file = paste("RFB_Accel_Models/Wind_selection_i/Wind_selection_interaction_gam_4k_", i, ".rds", sep = ""))
  
  }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# i = 10 looks marginally like it's the best
i = 10
full_mod_out <- readRDS(paste("RFB_Accel_Models/Wind_selection_i/Wind_selection_interaction_gam_4k_", i, ".rds", sep = ""))

# Subset data to include a certain number of points and weight the points by this value
# dat_outbound <- full_mod_out[["data"]]

# Check model out ####

pp_check(full_mod_out)

# Check that chains converge:
plot(full_mod_out)
# Check R-hat values:
summary(full_mod_out)

# Plot posterior intervals:
mcmc_plot(full_mod_out)

ce <- conditional_effects(full_mod_out)

ws_ws_ws_ce <- ce[["sc_ws"]]
ws_ce %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))

wd_ce <- ce[["sc_wd"]]
wd_ce %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))

wswd_ce <- ce[["sc_ws:sc_wd"]]
wswd_ce %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))

conditional_effects(full_mod_out)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT INTERACTION ####

# (Use original data)

new_plot_dat <- expand.grid(
  sc_ws = seq(min(dat_outbound$sc_ws, na.rm = TRUE), max(dat_outbound$sc_ws, na.rm = TRUE), length.out = 20),
  sc_wd = seq(min(dat_outbound$sc_wd, na.rm = TRUE), max(dat_outbound$sc_wd, na.rm = TRUE), length.out = 8)) %>%
  mutate(
    # Add predicted mean speed
    predicted_speed = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, mean),
    # Un-scale and un-center variables
    wind_speed = (sc_ws * sd(dat_outbound$wind_speed, na.rm = TRUE)) + mean(dat_outbound$wind_speed, na.rm = TRUE),
    relative_wind_dir = (sc_wd * sd(dat_outbound$relative_wind_dir, na.rm = TRUE)) + mean(dat_outbound$relative_wind_dir, na.rm = TRUE)) %>%
  mutate(
    # Add uncertainty intervals
    lower = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975))

# Summarise data to create histogram counts
h <- dat_outbound %>% group_by(Presence) %>%
  mutate(breaks = cut(wind_speed,
                      breaks = seq(0, 12, 0.5), 
                      labels = seq(0.25, 11.75, 0.5),
                      include.lowest = TRUE),
         breaks = as.numeric(as.character(breaks))) %>%
  group_by(Presence, breaks) %>% 
  summarise(n = n()) %>%
  mutate(pct = ifelse(Presence == 0, n/sum(n) * 0.6, 1 - n/sum(n) * 0.6)) 

dat_outbound <- dat_outbound %>%
  mutate(Presence_points = ifelse(Presence == 1, 1.01,
                                  ifelse(Presence == 0, -0.01, NA)))

mod_plot <- ggplot(dat_outbound) +
  # geom_segment(data = h, linewidth = 4, show.legend = FALSE,
  #              aes(x = breaks, xend = breaks, y = Presence, yend = pct, alpha = Presence)) +
  # scale_alpha_continuous(range = c(0.5,1)) +
  geom_point(aes(y = Presence_points, x = wind_speed, col = relative_wind_dir), alpha = 0.5) +
  geom_ribbon(data = new_plot_dat,
              aes(x = wind_speed, ymin = lower, ymax = upper, group = relative_wind_dir, fill = relative_wind_dir),
              alpha = 0.3) +
  geom_line(data = new_plot_dat, 
            aes(x = wind_speed, y = predicted_speed, group = relative_wind_dir, col = relative_wind_dir), linewidth = 1) +
  scale_colour_viridis_c(name = "Relative wind direction\n(tailwind \u2194 headwind)", breaks = c(0,45,90,135,180), labels = c(0,45,90,135,180), limits = c(0,180)) +
  scale_fill_viridis_c(name = "Relative wind direction\n(tailwind \u2194 headwind)", breaks = c(0,45,90,135,180), labels = c(0,45,90,135,180), limits = c(0,180)) +
  # facet_grid(Trip_section~Sex) +
  theme_bw() +
  labs(y = expression("Probability of selection"),
       x = expression("Wind speed (ms"^-1*")"),
       subtitle = "Outbound") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        legend.margin = margin(),
        axis.title = element_text(size = 8)) +
  guides(colour = guide_colorbar(barwidth = 15))

dat_outbound <- dat_outbound %>%
  mutate(Locations = ifelse(Presence == 1, "Original",
                            ifelse(Presence == 0, "Simulated", NA)))

speed_density <- ggplot(dat_outbound) +
  geom_density(aes(wind_speed, group = Locations, col = Locations, fill = Locations), alpha = 0.1) +
  scale_colour_manual(values = c("#0072B2", "grey")) +
  scale_fill_manual(values = c("#0072B2", "grey")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0,3,6,9,12), limits = c(0,12)) +
  labs(y = "Density",
       x = expression("Wind speed (ms"^-1*")"))

direction_density <- ggplot(dat_outbound) +
  geom_density(aes(relative_wind_dir, group = Locations, col = Locations, fill = Locations), alpha = 0.1) +
  scale_colour_manual(values = c("#0072B2", "grey")) +
  scale_fill_manual(values = c("#0072B2", "grey")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0,45,90,135,180), limits = c(0,180.1)) +
  scale_y_continuous(limits = c(0,0.008)) +
  labs(y = "Density",
       x = expression("Relative wind direction"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT WIND & ROTATED TACKS ####

# Get rotated map bits ready ####

chagos = read_sf("C:/Users/Ruth/Dropbox/HW_Seabird_Postdoc/RFB_Diving/RFB_Diving_Data/Chagos_Maps/chagos_maps/Chagos_v6_land_simple.shp")

orig_dat <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Wind_Calcs.csv") %>%
  filter(atCP == "No") %>%
  # Revalue 'Presence' column so that present = 1 and absent = 0
  mutate(Replicate = str_replace(as.character(replicate), "Replicate\\s?\\d+", "Replicate")) %>%
  mutate(Presence = as.numeric(dplyr::recode(Replicate, Replicate = "0", Original = "1"))) %>%
  na.omit() %>%
  filter(ID == "GY02207_trip2") %>%
  filter(Replicate == "Original") %>%
  mutate(Trip_section = ifelse(Trip_section == "Return", "Inbound", as.character(Trip_section))) %>%
  mutate(Trip_section = factor(Trip_section, levels = c("Outbound", "Middle", "Inbound")))

repl_dat <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Wind_Calcs.csv") %>%
  filter(atCP == "No") %>%
  # Revalue 'Presence' column so that present = 1 and absent = 0
  mutate(Replicate = str_replace(as.character(replicate), "Replicate\\s?\\d+", "Replicate")) %>%
  mutate(Presence = as.numeric(dplyr::recode(Replicate, Replicate = "0", Original = "1"))) %>%
  na.omit() %>%
  filter(ID == "GY02207_trip2") %>%
  mutate(Trip_section = factor(Trip_section, levels = c("Outbound", "Middle", "Return"))) %>%
  filter(Replicate == "Replicate") %>%
  filter(replicate == "Replicate 1" | replicate == "Replicate 2" | replicate == "Replicate 3" |
           replicate == "Replicate 4" | replicate == "Replicate 5" |
           replicate == "Replicate 6" | replicate == "Replicate 7" | replicate == "Replicate 8" |
           replicate == "Replicate 9" | replicate == "Replicate 10")

# Load background wind data ####

nc_dat <- nc_open("RFB_Accel_Data/Enviro_dat/Copernicus_2022_new/data_stream-oper.nc")

# Read in u and v wind speed data and wave height data (h)
# store the data in a 3-dimensional array
u_array <- ncvar_get(nc_dat, "u10") # 10 metre U wind component
v_array <- ncvar_get(nc_dat, "v10") # 10 metre V wind component

# Check what fill value is used for missing data:
fillvalue <- ncatt_get(nc_dat, "u10", "_FillValue")
# Replace these with NAs:
u_array[u_array == fillvalue$value] <- NA
v_array[v_array == fillvalue$value] <- NA

# Done, so close the netCDF files
nc_close(nc_dat)

rm(nc_dat, fillvalue)

# Create 3D raster bricks
u_brick <- brick(u_array,
                 xmn = min(repl_dat$Lon)-1, xmx = max(repl_dat$Lon)+1,
                 ymn = min(repl_dat$Lat)-1, ymx = max(repl_dat$Lat)+1,
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

v_brick <- brick(v_array,
                 xmn = min(repl_dat$Lon)-1, xmx = max(repl_dat$Lon)+1,
                 ymn = min(repl_dat$Lat)-1, ymx = max(repl_dat$Lat)+1,
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

u_mean <- calc(u_brick, fun = mean, na.rm = TRUE)
v_mean <- calc(v_brick, fun = mean, na.rm = TRUE)

rm(u_array, v_array, u_brick, v_brick)

# Create wind speed raster
wind_speed_raster <- sqrt(u_mean^2 + v_mean^2)

# Calculate wind direction:
wind_direction <- (180/pi) * atan2((u_mean/wind_speed_raster),(v_mean/wind_speed_raster))
wind_direction <- atan2(values(v_mean), values(u_mean)) * (180 / pi)
wind_direction <- (wind_direction + 360) %% 360  # Adjust to 0-360 degrees

# Create a new raster for wind direction:
wind_direction_raster <- wind_speed_raster  # Copy the structure of the zonal raster
values(wind_direction_raster) <- wind_direction

wind_df <- cbind(as.data.frame(wind_speed_raster, xy = TRUE),
                 as.data.frame(wind_direction_raster, xy = TRUE)[,3]) %>%
  rename(Lon = 1, Lat = 2, wind_speed = 3, wind_direction = 4) %>%
  mutate(u = cos(wind_direction * pi / 180),
         v = sin(wind_direction * pi / 180))

rm(wind_direction_raster, wind_speed_raster, wind_direction, u_mean, v_mean)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Show rotated tracks map ####

map_plot_a <- ggplot(repl_dat) +
  geom_raster(data = wind_df, aes(x = Lon, y = Lat, fill = wind_speed)) +
  geom_segment(data = wind_df, aes(x = Lon, y = Lat, xend = Lon + 0.2 * u, yend = Lat + 0.2 * v),
               arrow = arrow(length = unit(0.3, "cm")), linewidth = 0.3, alpha = 0.5) +
  scale_fill_viridis(option = "C", limits = c(6, 7), name = "Wind speed",
                     breaks = c(6,6.5,7), labels = c(6,6.5,7)) +
  # x = expression("Wind speed (ms"^-1*")")) +
  geom_path(data = orig_dat, aes(x = Lon, y = Lat), linewidth = 1) +
  geom_sf(data = chagos, fill = "#000000", col = "#000000") +
  theme_bw() +
  coord_sf(xlim = c(min(repl_dat$Lon), max(repl_dat$Lon)),
           ylim = c(min(repl_dat$Lat), max(repl_dat$Lat))) + # Restrict the visible area
  labs(x = "Lon", y = "Lat") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        legend.margin = margin(),
        axis.title = element_blank(),
        legend.key.height = unit(0.2, "cm")) +
  scale_x_continuous(breaks = pretty(repl_dat$Lon, n = 3)) +
  scale_y_continuous(breaks = pretty(repl_dat$Lat, n = 3))

map_plot_b <- ggplot(repl_dat) +
  geom_path(aes(x = Lon, y = Lat, group = "Replicate"),
            linewidth = 0.3, alpha = 0.5) +
  geom_path(data = orig_dat, aes(x = Lon, y = Lat,
                                 colour = Trip_section), linewidth = 1.2) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#CC79A7"), name = "") +
  geom_sf(data = chagos, fill = "#000000", col = "#000000") +
  theme_bw() +
  coord_sf(xlim = c(min(repl_dat$Lon), max(repl_dat$Lon)),
           ylim = c(min(repl_dat$Lat), max(repl_dat$Lat))) +
  labs(x = "Lon", y = "Lat") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        legend.margin = margin(),
        axis.title = element_blank()) +
  scale_x_continuous(breaks = pretty(repl_dat$Lon, n = 3)) +
  scale_y_continuous(breaks = pretty(repl_dat$Lat, n = 3))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What about the return leg? ####

return_mod <- readRDS("RFB_Accel_Models/Wind_selection_i/Wind_selection_interaction_return_gam_4k_10.rds")
summary(return_mod)

dat_inbound <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Wind_Calcs.csv") %>%
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
                            ifelse(Presence == 0, "Simulated", NA)))

new_plot_dat <- expand.grid(
  sc_ws = seq(min(dat_inbound$sc_ws, na.rm = TRUE), max(dat_inbound$sc_ws, na.rm = TRUE), length.out = 20),
  sc_wd = seq(min(dat_inbound$sc_wd, na.rm = TRUE), max(dat_inbound$sc_wd, na.rm = TRUE), length.out = 8)) %>%
  mutate(
    # Add predicted mean speed
    predicted_speed = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, mean),
    # Un-scale and un-center variables
    wind_speed = (sc_ws * sd(dat_inbound$wind_speed, na.rm = TRUE)) + mean(dat_inbound$wind_speed, na.rm = TRUE),
    relative_wind_dir = (sc_wd * sd(dat_inbound$relative_wind_dir, na.rm = TRUE)) + mean(dat_inbound$relative_wind_dir, na.rm = TRUE)) %>%
  mutate(
    # Add uncertainty intervals
    lower = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975))

# Summarise data to create histogram counts
h <- dat_inbound %>% group_by(Presence) %>%
  mutate(breaks = cut(wind_speed,
                      breaks = seq(0, 12, 0.5), 
                      labels = seq(0.25, 11.75, 0.5),
                      include.lowest = TRUE),
         breaks = as.numeric(as.character(breaks))) %>%
  group_by(Presence, breaks) %>% 
  summarise(n = n()) %>%
  mutate(pct = ifelse(Presence == 0, n/sum(n) * 0.6, 1 - n/sum(n) * 0.6)) 

return_mod_plot <- ggplot(dat_inbound) +
  # geom_segment(data = h, linewidth = 4, show.legend = FALSE,
  #              aes(x = breaks, xend = breaks, y = Presence, yend = pct, alpha = Presence)) +
  # scale_alpha_continuous(range = c(0.5,1)) +
  geom_point(aes(y = Presence_points, x = wind_speed, col = relative_wind_dir), alpha = 0.5) +
  geom_ribbon(data = new_plot_dat,
              aes(x = wind_speed, ymin = lower, ymax = upper, group = relative_wind_dir, fill = relative_wind_dir),
              alpha = 0.3) +
  geom_line(data = new_plot_dat, 
            aes(x = wind_speed, y = predicted_speed, group = relative_wind_dir, col = relative_wind_dir), linewidth = 1) +
  scale_colour_viridis_c(name = "Relative wind direction\n(tailwind \u2194 headwind)", breaks = c(0,45,90,135,180), labels = c(0,45,90,135,180), limits = c(0,180)) +
  scale_fill_viridis_c(name = "Relative wind direction\n(tailwind \u2194 headwind)", breaks = c(0,45,90,135,180), labels = c(0,45,90,135,180), limits = c(0,180)) +
  # facet_grid(Trip_section~Sex) +
  theme_bw() +
  labs(y = expression("Probability of selection"),
       x = expression("Wind speed (ms"^-1*")"),
       subtitle = "Inbound") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        legend.margin = margin(),
        axis.title = element_text(size = 8)) +
  guides(colour = guide_colorbar(barwidth = 15)) +
  scale_x_continuous(breaks = c(0,3,6,9,12), limits = c(0,12))

speed_density_in <- ggplot(dat_inbound) +
  geom_density(aes(wind_speed, group = Locations, col = Locations, fill = Locations), alpha = 0.1) +
  scale_colour_manual(values = c("#CC79A7", "grey")) +
  scale_fill_manual(values = c("#CC79A7", "grey")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0,3,6,9,12), limits = c(0,12)) +
  labs(y = "Density",
       x = expression("Wind speed (ms"^-1*")"))

direction_density_in <- ggplot(dat_inbound) +
  geom_density(aes(relative_wind_dir, group = Locations, col = Locations, fill = Locations), alpha = 0.1) +
  scale_colour_manual(values = c("#CC79A7", "grey")) +
  scale_fill_manual(values = c("#CC79A7", "grey")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0,45,90,135,180), limits = c(0,180.1)) +
  labs(y = "Density",
       x = expression("Relative wind direction"))

# return_mod_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Layout and save plots together ####

# design <- "AABB
#            AABB
#            CCCC
#            CCCC
#            CCCC
#            CCCC"
# 
# mod_plots <- mod_plot / return_mod_plot + plot_layout(guides = "collect", axes = "collect") &
#   theme(legend.position = "bottom")
# 
# map_plot_b + map_plot_a + mod_plots +
#   plot_layout(nrow = 2, design = design) +
#   plot_annotation(tag_levels = 'a') &
#   theme(legend.position = "bottom",
#         legend.box = "verticle",
#         legend.margin = margin(),
#         plot.tag = element_text(size = 10),
#         legend.title = element_text(size = 8))
# 
# ggsave("RFB_Accel_Plots/Wind_selectivity.png", height = 9, width = 7, dpi = 600)

# Try with new additions:

mod_plots <- mod_plot / return_mod_plot +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom")

density_plots <- speed_density / direction_density / speed_density_in / direction_density_in +
  plot_layout(axes = "collect") & theme(legend.position = "none")

design <- "AAAC
           BBBC"

selectivity <- mod_plots + density_plots + plot_layout(design = design)

(map_plot_b + map_plot_a) / selectivity +
  plot_layout(heights = c(1,2)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 10),
        legend.title = element_text(size = 8))

ggsave("RFB_Accel_Plots/Wind_selectivity_densities.png", height = 9, width = 7, dpi = 600)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Return trip bits for text:

ce <- conditional_effects(return_mod)

ws_ce <- ce[["sc_ws"]]
ws_ce %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))

wd_ce <- ce[["sc_wd"]]
wd_ce %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))

wswd_ce <- ce[["sc_ws:sc_wd"]]
wswd_ce %>%
  summarize(mean_effect = mean(estimate__),
            sd_effect = sd(estimate__),
            lower = mean(lower__),
            upper = mean(upper__))

AUC <- performance(prediction(predict(return_mod, type = "response")[,1],
                              as.vector(pull(na.omit(dati), Presence))),
                   measure = "auc")
AUC <- AUC@y.values[[1]]
