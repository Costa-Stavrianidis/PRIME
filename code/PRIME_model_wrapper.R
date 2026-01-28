# This wrapper script for running the model is adapted from the script
# provided for the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
#   https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/2e/Jags-Ybinom-XnomSsubjCcat-MbinomBetaOmegaKappa.R

source("utility_functions.R")

genMCMC = function(data, xName="x", NName="N", sName="s", gName="g", numSavedSteps=20000, 
                   saveName=NULL, thinSteps=1, runjagsMethod=runjagsMethodDefault, nChains=nChainsDefault) 
{ 
  require(rjags)
  require(runjags)
  
  # The data should have:
  # one variable x being a vector of integer # binomial successes (missense variants),
  # one variable N being a vector of integer # binomial trials (total variants),
  # one variable s being a factor of subregion identifiers,
  # and one variable g being a factor of gene identifiers, where
  # subregions are nested in genes.
  x = data[[xName]]
  N = data[[NName]]
  s = data[[sName]]
  g = data[[gName]]
  Nsubr = length(unique(s))
  Ngene =  length(unique(g))
  dataList = list(
    x = x,
    N = N,
    g = as.numeric(g),
    Nsubr = Nsubr,
    Ngene = Ngene
  )
  
  # Model code (this will create the PRIME_JAGS_model.txt file)
  modelString = "
  model {
    # Subregions
    for (sIdx in 1:Nsubr) {
        x[sIdx] ~ dbin(theta[sIdx], N[sIdx])
        theta[sIdx] ~ dbeta(omega[g[sIdx]] * (kappa[g[sIdx]] - 2) + 1, 
                            (1 - omega[g[sIdx]]) * (kappa[g[sIdx]] - 2) + 1)
    }

    # Genes
    for (gIdx in 1:Ngene) {
        omega[gIdx] ~ dbeta(omegaO * (kappaO - 2) + 1, 
                            (1 - omegaO) * (kappaO - 2) + 1)
        kappa[gIdx] <- kappaMinusTwo[gIdx] + 2
        kappaMinusTwo[gIdx] ~ dgamma(3, 0.1) T(0.001, )  # mean 30, truncated to avoid instability
    }

    # Genome
    omegaO ~ dbeta(1.0, 1.0)
    kappaO <- 10  # fixed concentration for genome prior
  }
  "
  writeLines(modelString , con = "PRIME_JAGS_model.txt")

  # Initialize chains
  initsList = function() {
    thetaInit = rep(NA, Nsubr)
    for (sIdx in 1:Nsubr) {  # for each subregion, resample data and compute proportion of missense
      resampledX = rbinom(1, size=N[sIdx], prob=x[sIdx]/N[sIdx])
      thetaInit[sIdx] = resampledX / N[sIdx]
    }
    thetaInit = 0.001 + 0.998 * thetaInit  # constrain to not be 0 or 1
    return(list(theta = thetaInit, 
                omega = aggregate(thetaInit, by=list(g), FUN=mean)$x,
                omegaO = mean(thetaInit)))
  }

  # Fit model
  parameters = c("theta", "omega", "kappa", "omegaO") 
  adaptSteps = 500  # Number of steps to adapt the samplers
  burnInSteps = 500  # Number of steps to burn-in the chains
  
  useRunjags = TRUE
  if (useRunjags) {
    runJagsOut <- run.jags(method = runjagsMethod,
                           model = "PRIME_JAGS_model.txt", 
                           monitor = parameters, 
                           data = dataList,  
                           inits = initsList, 
                           n.chains = nChains,
                           adapt = adaptSteps,
                           burnin = burnInSteps, 
                           sample = ceiling(numSavedSteps/nChains),
                           thin = thinSteps,
                           summarise = FALSE,
                           plots = FALSE)
    codaSamples = as.mcmc.list(runJagsOut)
  } else {
    # Create, initialize, and adapt the model:
    jagsModel = jags.model("PRIME_JAGS_model.txt", data = dataList, inits = initsList, 
                           n.chains = nChains, n.adapt = adaptSteps)
    # Burn-in:
    cat("Burning in the MCMC chain...\n")
    update(jagsModel, n.iter = burnInSteps)
    # Save the MCMC chain
    cat("Sampling final MCMC chain...\n")
    codaSamples = coda.samples(jagsModel, variable.names = parameters, 
                               n.iter = ceiling(numSavedSteps * thinSteps / nChains), 
                               thin = thinSteps)
  }  

  if (!is.null(saveName)) {
    save(codaSamples, file = paste(saveName, "MCMC.Rdata", sep=""))
  }
  return(codaSamples)
}


smryMCMC = function(codaSamples, compVal = 0.5, rope = NULL, 
                    diffSVec = NULL, diffCVec = NULL, 
                    compValDiff = 0.0, ropeDiff = NULL, 
                    saveName = NULL) {
  mcmcMat = as.matrix(codaSamples, chains=TRUE)
  summaryInfo = NULL
  rowIdx = 0
  # omega:
  for (parName in grep("omega", colnames(mcmcMat), value = TRUE)) {
    summaryInfo = rbind(summaryInfo, 
                        summarizePost(mcmcMat[ ,parName],
                                      compVal = compVal, ROPE = rope))
    rowIdx = rowIdx + 1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # kappa:
  for (parName in grep("kappa", colnames(mcmcMat), value = TRUE)) {
    summaryInfo = rbind(summaryInfo, 
                        summarizePost(mcmcMat[ ,parName],
                                      compVal = compVal, ROPE = rope))
    rowIdx = rowIdx + 1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # theta:
  for (parName in grep("theta", colnames(mcmcMat), value = TRUE)) {
    summaryInfo = rbind(summaryInfo, 
                        summarizePost(mcmcMat[ ,parName],
                                      compVal = compVal, ROPE = rope))
    rowIdx = rowIdx + 1
    rownames(summaryInfo)[rowIdx] = parName
  }
  
  # save:
  if (!is.null(saveName)) {
    write.csv(summaryInfo, file = paste(saveName, "SummaryInfo.csv", sep=""))
  }
  show(summaryInfo)
  return(summaryInfo)
}
