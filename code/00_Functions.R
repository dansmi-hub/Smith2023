# Functions

library(tidyverse)
library(Hmsc)

# S1 Compute abundance and other metrics in a function and plot.

s1 = function(community_matrix) {
  S = rowSums(community_matrix > 0, na.rm = TRUE)
  P = colMeans(community_matrix > 0, na.rm = TRUE)
  LA = rowSums(community_matrix, na.rm = TRUE)
  SA = colMeans(community_matrix, na.rm = TRUE)/P
  
  par(mfrow = c(2,2))
  hist(S, main = "Species Richness (S)")
  hist(P, main = "Species Prevelance (P)")
  hist(log(LA), main = "Log Local Abundance (LA)")
  hist(log(SA), main = "Log Species Abundance (SA)")
  par(mfrow = c(1,1))
  
}

# -------------------------------------------------------------------------

make_pres_ab = function(matrix) {
  matrix = ifelse(matrix > 1 , 1, 0)
  matrix = ifelse(is.na(matrix) , 0, matrix)
  return(matrix)
}

# model_timer -------------------------------------------------------------

# Supply with a HMSC model object and should give a rough estimate of model runtime
model_timer = function(model, samples = 10, gamma = TRUE) {
  
  gamma = as.logical(gamma)
  
  start <- Sys.time()
  sampleMcmc(model,
             thin = 1,
             samples = samples,
             verbose = 1,
             updater = list(GammaEta = gamma))
  
  elapsed = difftime(Sys.time(), start, units = "hours")
  
  multiplier_single_chain = 1000/samples
  
  message(paste("Model runtime for", samples, "samples :", round(elapsed, digits = 2), "hours\n"))
  
  
  message(paste0("Typical runtime for 1000 samples: ", 
             round(elapsed*multiplier_single_chain, digits = 2),
             "hours\n"))
  
  message(paste0("Typical runtime for normal model (thin = 100, samples = 1000): ",
             round(100*(elapsed*multiplier_single_chain), digits = 2),
             "hours\n\n",
             "Expect some overhead (~10%) for parralel computation"))
}


# -------------------------------------------------------------------------
# this runs several HMSC models in loop. Highly recommended to run smaller models first then larger models
# This code allows small models to run first and continue the fitting on large models so preliminary analysis can be done
# or can check if model is converged without having to restart the fitting process manually.

model_run_iter = function(model,
                          save_dir, 
                          test = FALSE, 
                          gamma = TRUE, 
                          unix = FALSE, 
                          nChains = 4,
                          nParallel = 1) {
  
  m = model
  
  # Defaults
  model.directory = file.path(save_dir)
  
  # Test or not?
  test.run = test
  
  if (test.run) {
    thin = c(1,2,3)
    samples = 100
    transient = 0
    verbose = thin
    nChains = 1
    nParallel = 1
  } else {
    thin = c(1,10,100,1000)
    samples = 1000
    transient = 5*thin
    verbose = thin
    nChains = nChains
    nParallel = nParallel
  }
  
  for(thin in thin)
  
  {
    # timing
    start <- Sys.time()
    
    # custom transient for each iteration
    transient = 50 * thin
    
    # run model
    m = sampleMcmc(
      m,
      thin = thin,
      samples = samples,
      transient = transient,
      nChains = nChains,
      # initPar = "fixed effects",
      nParallel = nParallel,
      verbose = thin,
      updater=list(GammaEta=gamma),
      useSocket = !as.logical(unix)
    )
    
    # save model 
    filename = file.path(
      model.directory,
      paste0(
        "model_chains_",
        as.character(nChains),
        "_samples_",
        as.character(samples),
        "_thin_",
        as.character(thin),
        ".rds"
      )
    )
    
    # save here
    readr::write_rds(m, file = filename, compress = "xz")
    
    # elapsed time generated for model here
    elapsed = as.character(round(difftime(Sys.time(), start, units = "hours"), digits = 2))
    
    # also generate filename for runtime same dir
    filename_elapsed = file.path(
      model.directory,
      paste0(
        "model_chains_",
        as.character(nChains),
        "_samples_",
        as.character(samples),
        "_thin_",
        as.character(thin),
        "_",
        elapsed,
        "_hours.txt"
      )
    )
    
    # save here
    writeLines(elapsed, filename_elapsed)
    
  }
}



# Merge model formulas ----------------------------------------------------

merge.formula <- function(form1, form2, ...){
  
  # get character strings of the names for the responses 
  # (i.e. left hand sides, lhs)
  lhs1 <- deparse(form1[[2]])
  lhs2 <- deparse(form2[[2]])
  if(lhs1 != lhs2) stop('both formulas must have the same response')
  
  # get character strings of the right hand sides
  rhs1 <- strsplit(deparse(form1[[3]]), " \\+ ")[[1]]
  rhs2 <- strsplit(deparse(form2[[3]]), " \\+ ")[[1]]
  
  # create the merged rhs and lhs in character string form
  rhs <- c(rhs1, rhs2)
  lhs <- lhs1
  
  # put the two sides together with the amazing 
  # reformulate function
  out <- reformulate(rhs, lhs)
  
  # set the environment of the formula (i.e. where should
  # R look for variables when data aren't specified?)
  environment(out) <- parent.frame()
  
  return(out)
}


# PA Mutate ---------------------------------------------------------------

pa_mutate <-
  function(x) {
    if_else(x > 0, as.integer(1), as.integer(0))
  }

# Cross Validation Functions ----------------------------------------------

## Effective Samples
ess = function(m) {
  mpost = convertToCodaObject(m)
  par(mfrow = c(3, 2))
  ess.beta = effectiveSize(mpost$Beta)
  psrf.beta = gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
  hist(ess.beta)
  hist(psrf.beta)
  ess.gamma = effectiveSize(mpost$Gamma)
  psrf.gamma = gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf
  hist(ess.gamma)
  hist(psrf.gamma)
  ns = m$ns
  sppairs = matrix(sample(x = 1:ns ^ 2, size = 100))
  tmp = mpost$Omega[[1]]
  for (chain in 1:length(tmp)) {
    tmp[[chain]] = tmp[[chain]][, sppairs]
  }
  ess.omega = effectiveSize(tmp)
  psrf.omega = gelman.diag(tmp, multivariate = FALSE)$psrf
  hist(ess.omega)
  hist(psrf.omega)
  # print("ess.rho:")
  # effectiveSize(mpost$Rho)
  # print("psrf.rho:")
  # gelman.diag(mpost$Rho)$psrf
}

## Fit Metrics
mfit = function(m) {
  preds = computePredictedValues(m)
  MF = evaluateModelFit(hM=m, predY=preds)
  par(mfrow = c(3,1))
  hist(MF$TjurR2, xlim = c(0,1), main=paste0("Mean TjurR2 = ", round(mean(MF$TjurR2),2)))
  hist(MF$AUC, xlim = c(0,1), main=paste0("Mean AUC = ", round(mean(MF$AUC),2)))
  hist(MF$RMSE, xlim = c(0,1), main=paste0("Mean RMSE = ", round(mean(MF$RMSE),2)))
  par(mfrow = c(1,1))
  return(MF)
}

## Cross Validated Fit Metrics
mfitcv = function(m, filepath) {
  output = list()
  # We can reduce thinning - accurate enough at 100 vs 1000 and saves 10x time
  m$thin = 100
  # Preds
  preds = computePredictedValues(m, verbose = m$thin)
  output$MF = evaluateModelFit(hM=m, predY=preds)
  # CV Preds
  folds = createPartition(m, 5)
  predsCV = computePredictedValues(m, partition = folds, verbose = m$thin)
  output$MFCV = evaluateModelFit(hM=m, predY=predsCV)
  pdf(paste0(filepath, ".pdf"), paper = "a4")
  par(mfrow = c(3,1))
  # Preds
  hist(output$MF$TjurR2, xlim = c(0,1), main=paste0("Mean TjurR2 = ", round(mean(output$MF$TjurR2),2)))
  hist(output$MF$AUC, xlim = c(0,1), main=paste0("Mean AUC = ", round(mean(output$MF$AUC),2)))
  hist(output$MF$RMSE, xlim = c(0,1), main=paste0("Mean RMSE = ", round(mean(output$MF$RMSE),2)))
  # CV
  hist(output$MFCV$TjurR2, xlim = c(0,1), main=paste0("Mean CV TjurR2 = ", round(mean(output$MFCV$TjurR2),2)))
  hist(output$MFCV$AUC, xlim = c(0,1), main=paste0("Mean CV AUC = ", round(mean(output$MFCV$AUC),2)))
  hist(output$MFCV$RMSE, xlim = c(0,1), main=paste0("Mean CV RMSE = ", round(mean(output$MFCV$RMSE),2)))
  par(mfrow = c(1,1))
  dev.off()
  
  return(output)
}







