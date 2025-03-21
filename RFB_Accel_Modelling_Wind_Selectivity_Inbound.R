rm(list=ls())

#### Habitat selection GLM of presence vs rotated data for the inbound trips ####

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
  na.omit()

ggplot(dat_inbound %>%
         mutate(year = lubridate::year(DateTime)) %>%
         filter(year == 2019)) +
  geom_histogram(aes(relative_wind_dir), bins = 15) +
  theme_bw() +
  labs(title = 2019) +
  ggplot(dat_inbound %>%
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
  dati <- dat_inbound %>%
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
  
  write_csv(as.data.frame(auc.vals), ("RFB_Accel_Models/Wind_selectivity_inbound_random_point_AUC_4k.csv"))

  # Save model:
  saveRDS(full_modi, file = paste("RFB_Accel_Models/Wind_selection_i/Wind_selection_inbound_interaction_gam_4k_", i, ".rds", sep = ""))
  
  }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# i = 10 looks like it's the best again
i = 10
full_mod_out <- readRDS(paste("RFB_Accel_Models/Wind_selection_i/Wind_selection_interaction_gam_4k_", i, ".rds", sep = ""))

# Subset data to include a certain number of points and weight the points by this value
# dat_inbound <- full_mod_out[["data"]]

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

