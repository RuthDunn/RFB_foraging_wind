#### Plot HMM results exploring the effects of wind and waves on RFB foraging ####

# Ruth Dunn 2025
# r.e.dunn@lancaster.ac.uk

rm(list=ls())

library(tidyverse)
library(data.table)
library(momentuHMM)
library(patchwork)
library(sf)
library(ggspatial) # to plot a scale bar

# Read in model and data:

model <- readRDS("RFB_Accel_Models/HMM_top_momentuHMM.rds")

gps_dat <- left_join(fread("RFB_Accel_Data/Accel_6_flight_variables/DG_GPS_segment_flight_info.csv"),
                     read_csv("RFB_Accel_Data/Accel_metadata.csv") %>%
                       mutate(Individual = ID) %>%
                       dplyr::select(Individual, Mass_deployment), by = "Individual") %>%
  mutate(group = with(rle(atCP), rep(seq_along(lengths), lengths))) %>%
  group_by(group) %>%
  filter(!(atCP == "Yes" & n() > 2)) %>%
  ungroup() %>%
  dplyr::select(-group) %>%
  mutate(wind_speed_sc = scale(wind_speed, center = T, scale = T),
         relwinddir_sc = scale(relative_wind_dir, center = T, scale = T),
         relwavedir_sc = scale(relative_wave_dir, center = T, scale = T),
         bodymass_sc = scale(Mass_deployment, center = T, scale = T))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot stationary state probabilities ####

stationary_ests <- plotStationary(model, plotCI = T, return = T, alpha = 0.95)

s_ests_wavedir <- stationary_ests$relwavedir_sc %>%
  imap_dfr(~ .x %>% mutate(States = .y)) %>%
  mutate(State = case_when(
    States == "state 1" ~ "Rest",
    States == "state 2" ~ "Forage",
    States == "state 3" ~ "Travel")) %>%
  mutate(var = cov * sd(gps_dat$relative_wave_dir, na.rm = T) +
           mean(gps_dat$relative_wave_dir, na.rm = T)) %>%
  mutate(Sig = "Non-significant")

s_ests_winddir <- stationary_ests$relwinddir_sc %>%
  imap_dfr(~ .x %>% mutate(States = .y)) %>%
  mutate(State = case_when(
    States == "state 1" ~ "Rest",
    States == "state 2" ~ "Forage",
    States == "state 3" ~ "Travel")) %>%
  mutate(var = cov * sd(gps_dat$relative_wind_dir, na.rm = T) +
           mean(gps_dat$relative_wind_dir, na.rm = T)) %>%
  mutate(Sig = "Non-significant")

s_ests_windspd <- stationary_ests$wind_speed_sc %>%
  imap_dfr(~ .x %>% mutate(States = .y)) %>%
  mutate(State = case_when(
    States == "state 1" ~ "Rest",
    States == "state 2" ~ "Forage",
    States == "state 3" ~ "Travel")) %>%
  mutate(var = cov * sd(gps_dat$wind_speed, na.rm = T) +
           mean(gps_dat$wind_speed, na.rm = T)) %>%
  mutate(Sig = ifelse(State == "Travel", "Non-significant", "Significant"))
    
s1 <- ggplot(s_ests_wavedir) +
  geom_line(aes(x = var, y = est, color = State, linetype = Sig), linewidth = 0.5, linetype = 2) +
  scale_colour_manual(values = c("#56B4E9", "#D55E00", "#F0E442")) +
  geom_ribbon(aes(x = var, ymin = lci, ymax = uci, fill = State), alpha = 0.3) +
  scale_fill_manual(values = c("#56B4E9", "#D55E00", "#F0E442")) +
  labs(y = "State occupancy probability",
       x = "Relative wave direction\n(receding waves \u2194 oncoming waves)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        title = element_text(size = 8)) +
  scale_x_continuous(limits = c(0, 180), breaks = c(0,45,90,135,180)) + guides(color = "none", fill = "none")

s2 <- ggplot(s_ests_winddir) +
  geom_line(aes(x = var, y = est, color = State), linewidth = 0.5, linetype = 2) +
  scale_colour_manual(values = c("#56B4E9", "#D55E00", "#F0E442")) +
  geom_ribbon(aes(x = var, ymin = lci, ymax = uci, fill = State), alpha = 0.3) +
  scale_fill_manual(values = c("#56B4E9", "#D55E00", "#F0E442")) +
  labs(y = "State occupancy probability",
       x = "Relative wind direction\n(tailwind \u2194 headwind)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        title = element_text(size = 8)) +
  scale_x_continuous(limits = c(0, 180), breaks = c(0,45,90,135,180)) + guides(color = "none", fill = "none")

s3 <- ggplot(s_ests_windspd) +
  geom_line(aes(x = var, y = est, color = State, linetype = Sig, linewidth = Sig)) +
  scale_colour_manual(values = c("#56B4E9", "#D55E00", "#F0E442"), guides) +
  geom_ribbon(aes(x = var, ymin = lci, ymax = uci, fill = State, alpha = Sig)) +
  scale_fill_manual(values = c("#56B4E9", "#D55E00", "#F0E442")) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_linewidth_manual(values = c(0.5, 1)) +
  scale_alpha_manual(values = c(0.3, 0.5)) +
  labs(y = "State occupancy probability",
       x = expression("Wind speed (ms"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 8),
        title = element_text(size = 8)) + guides(color = "none", fill = "none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# All tracks + HMM states ####

df_gps <- read.csv("RFB_Accel_Data/Accel_6_flight_variables/DG_HMMstates_TripSections.csv")
chagos = read_sf("RFB_Accel_Data/Chagos_Maps/Chagos_v6_land_simple.shp")

colonies <- tribble(
  ~Site,  ~C_Lat,    ~C_Lon,
  "East_Island",         -7.23,     72.42)

map <- ggplot(df_gps %>%
         mutate(State = ifelse(row_number() == 5246, "Transit", State)) %>%
         mutate(State = case_when(State == "Rest" ~ "Rest",
                                  State == "Forage" ~ "Feed",
                                  State == "Transit" ~ "Travel")) %>%
           mutate(State = factor(State, levels = c("Feed", "Travel", "Rest")))) +
  geom_path(aes(x = Lon, y = Lat, group = T, col = State)) +
  scale_colour_manual(values = c("#56B4E9", "#F0E442", "#D55E00")) +
  scale_x_continuous(limits = c(min(df_gps$Lon - 0.15), max(df_gps$Lon + 0.15)),
                     n.breaks = 5) +
  scale_y_continuous(limits = c(min(df_gps$Lat, na.rm = T), max(df_gps$Lat, na.rm = T))) +
  geom_sf(data = chagos, fill = "#000000", col = "#000000") +
  # geom_point(data = colonies, aes(x = C_Lon, y = C_Lat), col = "#FFB000", size = 3) +
  ylab("Latitude") +
  xlab("Longitude") +
  annotation_scale(location = "br", style = "ticks") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.justification = c(1, 0), legend.position = c(1-0.01, 0.05))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine the plots & save ####

a <- s3 + s2 + s1 +
  plot_layout(axes = "collect", guides = "collect") &
  theme(legend.position = "none")

map / a +
  plot_layout(axes = "collect", heights = c(2,1)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

ggsave("RFB_Accel_Plots/HMMs_stationary_probs.png", height = 8, width = 7, dpi = 600)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add icons ####

library(magick)

# Call back the plot
plot <- image_read("RFB_Accel_Plots/HMMs_stationary_probs.png")

# Import seabird images
RFB_travel <- image_read("RFB_Accel_Plots/booby_flying.png") %>%
  image_scale("400x400")
RFB_feed <- image_read("RFB_Accel_Plots/booby_diving.png") %>%
  image_scale("400x400")
RFB_rest <- image_read("RFB_Accel_Plots/booby_resting.png") %>%
  image_scale("400x400")

# Stack them on top of each other
plot1 <- image_composite(plot, RFB_travel, offset = "+3700+3650")
plot2 <- image_composite(plot1, RFB_feed, offset = "+3750+3200")
plot3 <- image_composite(plot2, RFB_rest, offset = "+3680+4100")

# And save
image_write(plot3, "RFB_Accel_Plots/Fig4_behavioural_states_RFBimage.png")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
