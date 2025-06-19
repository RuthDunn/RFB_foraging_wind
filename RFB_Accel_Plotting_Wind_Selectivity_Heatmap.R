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

dat_all <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Wind_Calcs.csv") %>%
  mutate(year = lubridate::year(DateTime)) %>%
  filter(replicate == "Original") %>%
  filter(Trip_section == "Outbound" | Trip_section == "Return") %>%
  group_by(ID) %>%
  # Extract the first, and last row of each trip:
  filter(row_number() == 1 | row_number() == n()) %>%
  ungroup()

summary(dat_all %>%
          filter(year == 2019) %>%
          dplyr::select(wind_direction))

summary(dat_all %>%
          filter(year == 2022) %>%
          dplyr::select(wind_direction))

ggplot(dat_all) +
  geom_histogram(aes(wind_direction), fill = "black", col = "black") +
  # geom_histogram((aes(bird_bearing, fill = Trip_section))) +
  # scale_fill_manual(values = c("#0072B2", "#CC79A7"), name = "") +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 360), limits = c(0,360), expand = c(0,0)) + # empty label for 360
  facet_grid(year ~ .) +
  coord_polar() +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())

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

# Get rotated map bits ready ####

chagos = read_sf("C:/Users/Ruth/Dropbox/HW_Seabird_Postdoc/RFB_Diving/RFB_Diving_Data/Chagos_Maps/chagos_maps/Chagos_v6_land_simple.shp")

# GY02207_trip2

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

repl_dat_out <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Wind_Calcs.csv") %>%
  filter(Trip_section == "Outbound" & atCP == "No") %>%
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

repl_dat_in <- read_csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections_Rotations_Inbound_Wind_Calcs.csv") %>%
  filter(Trip_section == "Return" & atCP == "No") %>%
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

nc_dat <- nc_open("RFB_Accel_Data/Enviro_dat/Copernicus_19_22_rotations/data_stream-oper_stepType-instant.nc")

# Read in times
t_wind <- ncvar_get(nc_dat, "valid_time")
t_dates <- c(seq(from = as.POSIXct("2019-02-03 00:00:00"), to = as.POSIXct("2019-02-19 23:00:00"), by = "hour"),
             seq(from = as.POSIXct("2022-02-03 00:00:00"), to = as.POSIXct("2022-02-19 23:00:00"), by = "hour"))

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
                 xmn = min(repl_dat_in$Lon)-1, xmx = max(repl_dat_out$Lon)+1,
                 ymn = min(repl_dat_out$Lat)-1, ymx = max(repl_dat_in$Lat)+1,
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

v_brick <- brick(v_array,
                 xmn = min(repl_dat_in$Lon)-1, xmx = max(repl_dat_out$Lon)+1,
                 ymn = min(repl_dat_out$Lat)-1, ymx = max(repl_dat_in$Lat)+1,
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Extract tracking period sections
times_df <- tibble(Layer = names(u_brick),
                   Index = t_wind,
                   DateTime = t_dates)

time_indices <- which(times_df$DateTime <= tail(orig_dat$DateTime,1) &
                        times_df$DateTime >= head(orig_dat$DateTime,1))
u_brick_time <- u_brick[[time_indices]]
v_brick_time <- v_brick[[time_indices]]

u_mean <- calc(u_brick_time, fun = mean, na.rm = TRUE)
v_mean <- calc(v_brick_time, fun = mean, na.rm = TRUE)

rm(u_array, v_array, u_brick, v_brick, u_brick_time, v_brick_time, time_indices, times_df)

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

# ROTATED TRACKS MAP ####

rotated_tracks_map <- ggplot(repl_dat_out) +
  geom_path(aes(x = Lon, y = Lat, group = "Replicate"),
            linewidth = 0.3, alpha = 0.5, col = "#0072B2") +
  geom_path(data = repl_dat_in,
            aes(x = Lon, y = Lat, group = "Replicate"),
            linewidth = 0.3, alpha = 0.5, col = "#CC79A7") +
  geom_path(data = orig_dat, aes(x = Lon, y = Lat,
                                 colour = Trip_section), linewidth = 1.2) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#CC79A7"), name = "") +
  geom_sf(data = chagos, fill = "#000000", col = "#000000") +
  theme_bw() +
  coord_sf(xlim = c(min(repl_dat_in$Lon, repl_dat_out$Lon), max(repl_dat_in$Lon, repl_dat_out$Lon)),
           ylim = c(min(repl_dat_in$Lat, repl_dat_out$Lat), max(repl_dat_in$Lat, repl_dat_out$Lat))) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        legend.margin = margin(),
        axis.title = element_blank()) +
  scale_x_continuous(breaks = pretty(repl_dat_in$Lon, n = 3)) +
  scale_y_continuous(breaks = pretty(repl_dat_in$Lat, n = 3)) +
  geom_point(aes(x = CP_Lon, y = CP_Lat), shape = 8, col = "#009E73")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# WIND EXAMPLE MAP ####

wind_example_map <- ggplot() +
  geom_raster(data = wind_df, aes(x = Lon, y = Lat, fill = wind_speed)) +
  geom_segment(data = wind_df, aes(x = Lon, y = Lat, xend = Lon + 0.2 * u, yend = Lat + 0.2 * v),
               arrow = arrow(length = unit(0.3, "cm")), linewidth = 0.3, alpha = 0.5) +
  scale_fill_viridis(option = "D", name = expression("Wind\nspeed (ms"^-1*")"))  +
  geom_path(data = orig_dat, aes(x = Lon, y = Lat), linewidth = 1) +
  geom_sf(data = chagos, fill = "#000000", col = "#000000") +
  theme_bw() +
  coord_sf(xlim = c(min(repl_dat_in$Lon, repl_dat_out$Lon), max(repl_dat_in$Lon, repl_dat_out$Lon)),
           ylim = c(min(repl_dat_in$Lat, repl_dat_out$Lat), max(repl_dat_in$Lat, repl_dat_out$Lat))) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        legend.margin = margin(),
        axis.title = element_blank(),
        legend.key.height = unit(0.2, "cm")) +
  scale_x_continuous(breaks = pretty(repl_dat_in$Lon, n = 3)) +
  scale_y_continuous(breaks = pretty(repl_dat_in$Lat, n = 3)) +
  geom_point(data = repl_dat_in, aes(x = CP_Lon, y = CP_Lat), shape = 8, col = "#009E73")

rotated_tracks_map + wind_example_map

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT INTERACTION (OUTBOUND) ####

full_mod_out <- readRDS(paste("RFB_Accel_Models/Wind_selection_i/Wind_selection_interaction_gam_4k_", 10, ".rds", sep = ""))

# Create grid in original units, then manually scale
new_plot_dat <- expand.grid(
  wind_speed = seq(0, 12, length.out = 12),
  relative_wind_dir = seq(0, 180, length.out = 14)
) %>%
  # Manually scale using the same parameters as in your model
  mutate(
    sc_ws = (wind_speed - mean(dat_outbound$wind_speed, na.rm = TRUE)) / sd(dat_outbound$wind_speed, na.rm = TRUE),
    sc_wd = (relative_wind_dir - mean(dat_outbound$relative_wind_dir, na.rm = TRUE)) / sd(dat_outbound$relative_wind_dir, na.rm = TRUE)
  ) %>%
  # Predict from the model
  mutate(
    predicted_speed = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, mean),
    lower = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(full_mod_out, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975)
  ) %>%
  # Bin wind speeds for plotting
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 12, by = 1),
                              labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6",
                                         "6-7", "7-8", "8-9", "9-10", "10-11", "11-12"),
                              include.lowest = TRUE,
                              right = FALSE))

dat_outbound <- dat_outbound %>%
  mutate(Presence_points = ifelse(Presence == 1, 1.01,
                                  ifelse(Presence == 0, -0.01, NA)))
mod_plot <- ggplot(new_plot_dat, 
                       aes(x = relative_wind_dir, 
                           y = wind_speed_bin, 
                           fill = predicted_speed)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", name = "Probability of selection  ",
                       limits = c(0,1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  labs(x = expression("Relative wind direction (tailwind \u2194 headwind)"),
       y = expression("Wind speed (ms"^-1*")"),
       subtitle = "Outbound") +
  theme(axis.title = element_text(size = 8),
        legend.position = "right") +
  scale_x_continuous(breaks = c(0, 45, 90, 135, 180))

mod_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# OUTBOUND WIND DENSITY PLOTS ####

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
       x = expression("Wind speed"))

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

direction_density/speed_density

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
  mutate(Presence_points = ifelse(Presence == 1, 1.01,
                                  ifelse(Presence == 0, -0.01, NA))) %>%
  mutate(Locations = ifelse(Presence == 1, "Original",
                            ifelse(Presence == 0, "Simulated", NA))) %>%
  mutate(wind_speed_bin = cut(wind_speed, 
                              breaks = seq(0, 11, by = 1),  # Defines bin edges
                              labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6",
                                         "6-7", "7-8", "8-9", "9-10", "10-11"),  # Labels for bins                              include.lowest = TRUE,  # Ensures 0 is included
                              right = FALSE))

# Use original scale (not scaled wind direction)
new_plot_dat_in <- expand.grid(
  wind_speed = seq(0, 11, length.out = 11),
  relative_wind_dir = seq(0, 180, length.out = 14)
) %>%
  # Now scale manually using the same parameters as in the model
  mutate(
    sc_ws = (wind_speed - mean(dat_inbound$wind_speed, na.rm = TRUE)) / sd(dat_inbound$wind_speed, na.rm = TRUE),
    sc_wd = (relative_wind_dir - mean(dat_inbound$relative_wind_dir, na.rm = TRUE)) / sd(dat_inbound$relative_wind_dir, na.rm = TRUE)
  ) %>%
  mutate(
    predicted_speed = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, mean),
    lower = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.025),
    upper = posterior_epred(return_mod, newdata = ., re_formula = NA) %>% apply(2, quantile, probs = 0.975),
    wind_speed_bin = cut(wind_speed,
                         breaks = seq(0, 12, by = 1),
                         labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6",
                                    "6-7", "7-8", "8-9", "9-10", "10-11", "11-12"),
                         right = FALSE)
  )

dat_inbound <- dat_inbound %>%
  mutate(Presence_points = ifelse(Presence == 1, 1.01,
                                  ifelse(Presence == 0, -0.01, NA)))

return_mod_plot <- ggplot(new_plot_dat_in, 
                          aes(x = relative_wind_dir, 
                              y = wind_speed_bin, 
                              fill = predicted_speed)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", name = "Probability of selection  ",
                       limits = c(0,1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  labs(x = expression("Relative wind direction (tailwind \u2194 headwind)"),
       y = expression("Wind speed (ms"^-1*")"),
       subtitle = "Inbound") +
  theme(axis.title = element_text(size = 8),
        legend.position = "right") +
  scale_x_continuous(breaks = c(0, 45, 90, 135, 180))

return_mod_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# INBOUND WIND DENSITY PLOTS ####

speed_density_in <- ggplot(dat_inbound) +
  geom_density(aes(wind_speed, group = Locations, col = Locations, fill = Locations), alpha = 0.1) +
  scale_colour_manual(values = c("#CC79A7", "grey")) +
  scale_fill_manual(values = c("#CC79A7", "grey")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0,3,6,9,12), limits = c(0,12)) +
  labs(y = "Density",
       x = expression("Wind speed"))

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

mod_plots <- mod_plot / return_mod_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  guides (fill = guide_colourbar(barwidth = 10, barheight = 0.5))

density_plots <- direction_density / speed_density / direction_density_in / speed_density_in +
  plot_layout(axes = "collect") & theme(legend.position = "none")

design <- "AAAC
           BBBC"

selectivity <- mod_plots + density_plots + plot_layout(design = design)

(rotated_tracks_map + wind_example_map) / selectivity +
  plot_layout(heights = c(1,2)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        title = element_text(size = 8))

ggsave("RFB_Accel_Plots/Wind_selectivity_heatmap_new.png", height = 9, width = 7, dpi = 600)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add booby image to plot ####

library(magick)

# Call back the plot
plot <- image_read("RFB_Accel_Plots/Wind_selectivity_heatmap_new.png")

# Import seabird images
RFB_img <- image_read("RFB_Accel_Plots/booby_flying.png") %>%
  image_scale("700x700")

# Stack them on top of each other
plot1 <- image_composite(plot, RFB_img, offset = "+3400+60")

# And save
image_write(plot1, "RFB_Accel_Plots/Final/Fig2_heatmap_RFBimage.png")

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

