# These utility functions for running the model are adapted from the utility functions
# provided for the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
#   https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/2e/DBDA2E-utilities.R

# Ensures correct packages are installed
want = c("parallel", "rjags", "runjags", "compute.es")
have = want %in% rownames(installed.packages())
if (any(!have)) {install.packages(want[!have])}

# Load necessary packages
try(library(rjags))
try(library(runjags))
try(runjags.options(inits.warning = FALSE, rng.warning = FALSE))

# Running in parallel
library(parallel)
nCores = detectCores() 
if (!is.finite(nCores)) {nCores = 1} 
# Parallel options
if (nCores > 4) { 
  nChainsDefault = 4
  runjagsMethodDefault = "parallel"
}
if (nCores == 4) { 
  nChainsDefault = 3
  runjagsMethodDefault = "parallel"
}
# Non-parallel
if (nCores < 4) { 
  nChainsDefault = 3 
  runjagsMethodDefault = "rjags"
}

# HDI

HDIofMCMC = function(sampleVec, credMass = 0.95) {
  sortedPts = sort(sampleVec)
  ciIdxInc = ceiling(credMass * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
}

# Summary statistics

summarizePost = function(paramSampleVec, 
                        compVal = NULL, ROPE = NULL, credMass = 0.95) {
  meanParam = mean(paramSampleVec)
  medianParam = median(paramSampleVec)
  dres = density(paramSampleVec)
  modeParam = dres$x[which.max(dres$y)]
  mcmcEffSz = round(effectiveSize(paramSampleVec), 1)
  names(mcmcEffSz) = NULL
  hdiLim = HDIofMCMC(paramSampleVec, credMass = credMass)
  if (!is.null(compVal)) {
    pcgtCompVal = (100 * sum(paramSampleVec > compVal) 
                    / length(paramSampleVec))
  } else {
    compVal = NA
    pcgtCompVal = NA
  }
  if (!is.null(ROPE)) {
    pcltRope = (100 * sum(paramSampleVec < ROPE[1]) 
                 / length(paramSampleVec))
    pcgtRope = ( 100 * sum(paramSampleVec > ROPE[2]) 
                 / length(paramSampleVec))
    pcinRope = 100 - (pcltRope + pcgtRope)
  } else { 
    ROPE = c(NA, NA)
    pcltRope = NA 
    pcgtRope = NA 
    pcinRope = NA 
  }  
  return(c(Mean = meanParam, Median = medianParam, Mode = modeParam, 
             ESS = mcmcEffSz,
             HDImass = credMass, HDIlow = hdiLim[1], HDIhigh = hdiLim[2], 
             CompVal = compVal, PcntGtCompVal = pcgtCompVal, 
             ROPElow = ROPE[1], ROPEhigh = ROPE[2],
             PcntLtROPE = pcltRope, PcntInROPE = pcinRope, PcntGtROPE = pcgtRope))
}