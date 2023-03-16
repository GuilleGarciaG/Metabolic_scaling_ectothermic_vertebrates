#simulate one data set under model with student-t errors
#Arguments:
#realdata: data frame of real data, containing at least the variables group (water-breather or air-breather), log10_L, experiment
#trueparms: list of true parameters from which to simulate: 
#b_Intercept, b_log10_L, b_groupwaterMbreather, b_log10L:groupwaterMbreather, sigma (scale for t errors), mu (location for t errors)
#sd_experiment__Intercept and sd_experiment_log10_L: standard deviations for random effects of experiment on slope and intercept
#cor_experiment__Intercept__log10_L: correlation between slope and intercept experiment-level effects
#sd_phylo__Intercept: standard deviation for phylogenetic effects on intercept
#A: known phylogenetic covariance matrix
#Value: a simulated data set with same size and structure as real data set
simulatetdata <- function(realdata, trueparms, A){
  simdata <- realdata
  simdata$group <- as.numeric(simdata$group == "water-breather")
  nexperiment <- length(unique(simdata$experiment))
  excov <- matrix(c(trueparms$sd_experiment__Intercept^2, trueparms$cor_experiment__Intercept__log10_L * trueparms$sd_experiment__Intercept * trueparms$sd_experiment__log10_L, trueparms$cor_experiment__Intercept__log10_L * trueparms$sd_experiment__Intercept * trueparms$sd_experiment__log10_L, trueparms$sd_experiment__log10_L^2), nrow = 2, ncol = 2) #covariance matrix for experiment effects (slope and intercept drawn from bivariate normal)
  
  eeffects <- rmvnorm(n = nexperiment, mean = c(0, 0), sigma = excov) #one draw per experiment: first column is intercept, second column  is slope
  eindex <- as.numeric(factor(simdata$experiment)) #index the experiments
  simdata$ei <- eeffects[eindex, 1] #select the experiment intercepts
  simdata$es <- eeffects[eindex, 2] #select the experiment slopes
  
  nphylo <- length(unique(simdata$phylo))
  phyloeffect <- rmvnorm(n = 1, mean = rep(0, nphylo), sigma = trueparms$sd_phylo__Intercept^2 * cov2cor(A)) #one draw from multivariate normal witrh covariance matrix proportional to A
  pindex <- as.numeric(factor(simdata$phylo))
  simdata$ep <- phyloeffect[pindex] #select phylogenetic effects 
  mu <- trueparms$b_Intercept + trueparms$b_log10_L * simdata$log10_L + trueparms$b_groupwaterMbreather * simdata$group + trueparms$`b_log10_L:groupwaterMbreather` * simdata$log10_L * simdata$group + simdata$ei + simdata$es * simdata$log10_L + simdata$ep

  simdata$b <- rt(n = dim(simdata)[1], df = trueparms$nu) * trueparms$sigma + mu #sample from non-standardized t-distribution
  return(simdata)
}

#read files containing models fitted to simulated data and plot the posterior distributions
#Arguments:
#path: path to directory containing fitted models
#fnames: vector of file names of fitted models
#trueparms: list of true parameter values
#Value: for each parameter, plot posterior density for each simulated data set, and indicate true value
plotsimulated <- function(path, fnames, trueparms){
  niter <- length(fnames)
  nparms <- length(trueparms)
  par(mfrow = c(3, 4))
  for(i in 1:nparms){
    whichparm <- names(trueparms)[i]
    simfit <- readRDS(file = paste(path, fnames[1], sep = ""))
    estimated <- as_draws_matrix(simfit, variable = whichparm)
    plot(density(estimated), col = adjustcolor("black", 0.4), xlab = whichparm, main = "")
    for(j in 2:niter){
      simfit <- readRDS(file = paste(path, fnames[j], sep = ""))
      estimated <- as_draws_matrix(simfit, variable = whichparm)
      lines(density(estimated), col = adjustcolor("black", 0.2))
    }
    abline(v = trueparms[i], lwd = 2)
  }
}