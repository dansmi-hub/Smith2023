## ---------------------------

library(Hmsc)

# -------------------------------------------------------------------------

source("code/00_Functions.R")

# model_run ---------------------------------------------------------------

model_run_iter(model = read_rds("data/Model.rds"), 
               save_dir = "results/fit/", 
               test = FALSE, unix = TRUE, nChains = 4, nParallel = 4)
