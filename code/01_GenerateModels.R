## 01_Generate_Models


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(Hmsc)

# Load --------------------------------------------------------------------

data <- read_rds("data/prepared_data.rds")

# Generate unique XY for each plot_ID
xy <- distinct(data.frame(data$xy))
rownames(xy) <- unique(data$plot_id)

# Studydesign
StudyDesign <- data.frame(
  plot_id = as.factor(data$plot_id),
  year = as.factor(paste0("Year_", data$year)),
  season = as.factor(paste0("Season_", data$season))
)

# RL
rL.year <- HmscRandomLevel(units = unique(StudyDesign$year))
rL.season <- HmscRandomLevel(units = unique(StudyDesign$season))
rL.spatial <- HmscRandomLevel(sData = xy)

# All random levels
rl.full <- list(
  plot_id = rL.spatial,
  season = rL.season,
  year = rL.year
)

# X formula
xform = ~ width +
  water_temp +
  turbidity +
  salinity +
  ph +
  height_vertical +
  height_emergent +
  cover_floating +
  cover_vertical +
  cover_emergent +
  shaded +
  dO2

model = Hmsc(
  # This is a full model with co-occurring predators
  # A matrix of species occurrence/abundance records
  Y = make_pres_ab(as.matrix(data$Y)),
  # A data frame of measured covariates
  XData = data$X,
  # A formula object for linear regression
  XFormula = xform,
  # Distribution
  distr = "probit",
  # Study Design
  studyDesign = StudyDesign,
  # Random Levels
  ranLevels = rl.full
)

# Save The New Models -------------------------------------------------------

write_rds(model, "data/Model.rds")

# Run Model -----------------------------------------------------------------

model_run_iter(model = read_rds("data/Model.rds"), 
               save_dir = "results/fit/", 
               test = FALSE, unix = TRUE, nChains = 4, nParallel = 4)

# Cross Validate Models -----------------------------------------------------

modelcv = mfitcv(model)
write_rds(modelcv, "results/fit/Model_CV.rds")














