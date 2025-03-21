#### Run HMMs to explore the effect of wind and waves on RFB foraging ####

# Ruth Dunn 2025
# r.e.dunn@lancaster.ac.uk

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get ready ####

# install.packages("rngtools")
# install.packages("momentuHMM")

rm(list=ls())

# Set CRS:
original_crs <- 4326
new_crs <- 32643

library(tidyverse)
library(data.table)
library(sf) # for spatial stuff
library(momentuHMM) # to run HMMs over GPS data
library(slider) # for rolling calculations
library(patchwork) # for combining plots
library(ncdf4) # package for netcdf manipulation
library(raster) # to work with raster data
library(patchwork) # for saving plots

# Read in the interpolated data that has wind and wave data associated with it

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

# HMM Bit: ####

# Prep data for moveHMM

hmmdata_prep <- data.frame(x = gps_dat$Lon,
                           y = gps_dat$Lat,
                           ID = gps_dat$ID,
                           wind_speed_sc = gps_dat$wind_speed_sc,
                           relwinddir_sc = gps_dat$relwinddir_sc,
                           relwavedir_sc = gps_dat$relwavedir_sc,
                           bodymass_sc = gps_dat$bodymass_sc)

hmmdata <- prepData(hmmdata_prep, type = "LL", coordNames = c("x", "y"),
                    covNames=c("wind_speed_sc", "relwinddir_sc", "relwavedir_sc", "bodymass_sc"))

# Use step lengths and turning angles from moveHMM model:

# old_mod <- readRDS("RFB_Accel_Models/HMM_moveHMM.rds")

stepMean0 <- c(0.10315611, 0.7220658, 3.182635)
stepSD0 <- c(0.04777501, 0.7192329, 0.640110)

angleMean0 <- c(-0.003215519, -0.02672849, -0.0009539601)
angleCon0 <- c(23.159656351, 0.82547652, 19.3683480639)

stepPar0 <-c(stepMean0, stepSD0)
anglePar0 <-c(angleMean0, angleCon0)

# Transition probability formula:
formula <- list()

formula[[2]] <- ~ wind_speed_sc : relwinddir_sc + wind_speed_sc : relwavedir_sc + wind_speed_sc + relwinddir_sc + relwavedir_sc

formula[[3]] <- ~ wind_speed_sc : relwinddir_sc + relwavedir_sc + wind_speed_sc + relwinddir_sc
formula[[4]] <- ~ wind_speed_sc : relwavedir_sc + relwinddir_sc + wind_speed_sc + relwavedir_sc

formula[[5]] <- ~ wind_speed_sc + relwinddir_sc + relwavedir_sc

formula[[6]] <- ~ wind_speed_sc + relwinddir_sc
formula[[7]] <- ~ wind_speed_sc + relwavedir_sc

formula[[8]] <- ~ wind_speed_sc
formula[[9]] <- ~ relwinddir_sc
formula[[10]] <- ~ relwavedir_sc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run through model combinations:

m.list <- vector(mode = "list", length =length(formula))
for (i in 2:length(formula)) {
  print(i)
  m.list[[i]] <- fitHMM(data = hmmdata, nbStates = 3,
                        dist = list(step = "gamma", angle = "vm"),
                        Par0 = list(step = stepPar0, angle = anglePar0),
                        estAngleMean = list(angle=TRUE),
                        formula = formula[[i]])
  model <- m.list[[i]]
}

# Check which model is best:

out.df <- vector(mode = "list", length =length(formula))
for (i in 2:length(formula)) {
  print(i)
  file.in <- paste0("./RFB_Accel_Models/", paste0("HMM_mod_", i, ".RData"))
  if (i == 1) { model <- m1} else { model <- model}
  # outputting AICs into dataframe
  if (i == 1) { form_out <- 1} else { form_out <- as.character(formula[[i]])[2]}
  load(file = file.in)
  out.df[[i]] <- data.frame(Model = paste0("M", i),
                            Formula = form_out, AIC = AIC(model))}

all.out <- do.call(rbind, out.df)
all.out <- all.out[order(all.out$AIC),]

all.out

all.out$AIC[2] - all.out$AIC[1]

write.csv(all.out, "RFB_Accel_Models/HMM_AIC_comparisons.csv")
all.out <- read_csv("RFB_Accel_Models/HMM_AIC_comparisons.csv")

# 1 M2 wind_speed_sc:relwinddir_sc + wind_speed_sc:relwavedir_sc + wind_speed_sc + relwinddir_sc + relwavedir_sc 30698.56

i = 1 + 1
model <- m.list[[i]]

saveRDS(m.list[[i]], "RFB_Accel_Models/HMM_top_momentuHMM.rds")
model <- readRDS("RFB_Accel_Models/HMM_top_momentuHMM.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check the model out ####

plot(model, plotCI = TRUE)

# State 1 (yellow): small step lengths (slow), narrow turning angles = REST
# State 2 (blue): small step lengths (slow), wide turning angles = FORAGE
# State 3 (green): large step lengths (fast), narrow turning angles = TRANSIT

# Decode most likely state sequence
states <- viterbi(model)
# Derive percentage of time spent in each state
table(states)/nrow(gps_dat)

# Create transitions by pairing each element with the next
transitions <- data.frame(
  from = head(states, -1),
  to = tail(states, -1)
)
# Count occurrences of each transition
table(transitions)

# Access specific transitions in the table
nrow(subset(transitions, from == 2 & to == 1)) # F -> R
nrow(subset(transitions, from == 2 & to == 3)) # F -> T
nrow(subset(transitions, from == 1 & to == 2)) # R -> F
nrow(subset(transitions, from == 3 & to == 2)) # T -> F

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pull out values for text ####

beta.full <- CIbeta(model)$beta
beta.full.est <- as.data.frame(beta.full$est)
beta.full.upr <- as.data.frame(beta.full$upper)
beta.full.lwr <- as.data.frame(beta.full$lower)
beta.df <- data.frame(Est = c(beta.full.est$`1 -> 2`, beta.full.est$`1 -> 3`,
                              beta.full.est$`2 -> 1`, beta.full.est$`2 -> 3`,
                              beta.full.est$`3 -> 1`, beta.full.est$`3 -> 2`), 
                      Upr = c(beta.full.upr$`1 -> 2`, beta.full.upr$`1 -> 3`,
                              beta.full.upr$`2 -> 1`, beta.full.upr$`2 -> 3`,
                              beta.full.upr$`3 -> 1`, beta.full.upr$`3 -> 2`), 
                      Lwr = c(beta.full.lwr$`1 -> 2`, beta.full.lwr$`1 -> 3`,
                              beta.full.lwr$`2 -> 1`, beta.full.lwr$`2 -> 3`,
                              beta.full.lwr$`3 -> 1`, beta.full.lwr$`3 -> 2`), 
                      Transitions = rep(colnames(beta.full.est), each = nrow(beta.full.est)),
                      Covariates = rep(rownames(beta.full.est), 3)) %>%
  filter(Covariates != "(Intercept)") %>%
  filter(Transitions != "1 -> 3" &
           Transitions != "3 -> 1") %>%
  mutate(Transitions = case_when(Transitions == "1 -> 2" ~ "R -> F",
                                 Transitions == "2 -> 3" ~ "F -> T",
                                 Transitions == "2 -> 1" ~ "F -> R",
                                 Transitions == "3 -> 2" ~ "T -> F")) %>%
  mutate(Covariates = case_when(Covariates == "relwinddir_sc" ~ "Relative wind direction",
                                Covariates == "relwavedir_sc" ~ "Relative wave direction",
                                Covariates == "wind_speed_sc" ~ "Wind speed")) %>%
  mutate(Covariates = factor(Covariates, levels = desired_order))

# Wind speed on transit to forage

beta.df %>%
  filter(Covariates == "Wind speed") %>%
  filter(Transitions == "T -> F")

beta.df %>%
  filter(Covariates == "Wind speed") %>%
  filter(Transitions == "F -> T")

# Wind direction on transit to forage

beta.df %>%
  filter(Covariates == "Relative wind direction") %>%
  filter(Transitions == "T -> F")

beta.df %>%
  filter(Covariates == "Relative wind direction") %>%
  filter(Transitions == "F -> T")

# Wave direction on transit to forage

beta.df %>%
  filter(Covariates == "Relative wave direction") %>%
  filter(Transitions == "T -> F")

beta.df %>%
  filter(Covariates == "Relative wave direction") %>%
  filter(Transitions == "F -> T")

# Forage to rest in tailwinds

beta.df %>%
  filter(Covariates == "Relative wind direction") %>%
  filter(Transitions == "F -> R")

# Forage to rest in low wind speeds

beta.df %>%
  filter(Covariates == "Wind speed") %>%
  filter(Transitions == "F -> R")

# Forage to rest in low wind speeds

beta.df %>%
  filter(Covariates == "Relative wind direction") %>%
  filter(Transitions == "R -> F")

