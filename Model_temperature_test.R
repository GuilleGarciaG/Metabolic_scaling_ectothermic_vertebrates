################################################################################
# Script to fit phylogenetic multilevel model using brms: Part 1
# Title: Temperature effect on metabolic scaling of ectothermic vertebrates
# Authors: Guillermo García-Gómez (guillegar.gz@gmail.com) & Matthew Spencer
# Date: 16032023
# Operating System: Windows 10 Pro 21H2
# ------------------------------------------------------------------------------
# Cite as: ?
# ------------------------------------------------------------------------------
rm(list=ls())# clear the work environment
today <- format(Sys.Date(),"%Y%m%d")# setting the date
# ------------------------------------------------------------------------------
# It needs to be set to Project directory
getwd()# to check
#
# ------------------------------------------------------------------------------
# Packages ####
#
# Bayesian analysis
library(rstan)
library(rstantools)
library(brms)
library(bayestestR)
library(loo)

# visual-checking model performance
library(ggmcmc)
library(bayesplot)
library(tidybayes)

# phylogenetic tree
library(rotl)
library(ape)

## plotting 
library(ggplot2)
library(ggtree)
library(cowplot)

# data sorting
library(dplyr)

# Computes t probabilities, quantiles, random deviates and densities
library(mvtnorm)

# call functions to check model performance
source("temperaturefuncs.R")

fitt <- TRUE # fit model with student-t errors?
fitgaussian <- FALSE # fit model with Gaussian errors?
fittsimulated <- FALSE # fit student-t model to data simulated under posterior mean parameters estimated from real data?
tdiagnostics <- FALSE # diagnostics for model with student-t errors
processtsimulated <- FALSE # plot estimates from student-t model fitted to simulated data
simfitname <- "../Model_outputs/simfitt" # base for names in which saved fits to simulated data will be stored

# ------------------------------------------------------------------------------
set.seed(1234)# to replicate the results

# General STAN specifications
rstan::rstan_options(auto_write = TRUE) # write objects to the hard disk to avoid recompilation 
options(mc.cores = parallel::detectCores())

#### Load data ####
all_temp <- read.csv("table_S1.csv") # Table S1

# Histograms: 

# response variable: slopes b
hist(all_temp$b)

# predictor variable: log10(Metabolic level, L)
hist(log10(all_temp$L))

#### Models for all species (water + air-breathers) ####

# Building the phylogenetic tree for all species (needs internet connection)

all_temp.taxa <- tnrs_match_names(as.character(all_temp$species_phylo)) # match species from dataset to those included in Tree of Life
all_temp.tree <- tol_induced_subtree(ott_ids = ott_id(all_temp.taxa), label = "name") # building a tree

# Check phylogenetic tree
ggtree(all_temp.tree, layout = "fan", open.angle=0, ladderize = TRUE, size=0.75) +
  geom_tiplab(size=2.1, hjust = -.1, color='blue')

# group species colour-coded according to respiration mode by selecting nodes in the tree 
(tree2 <- groupClade(all_temp.tree, c(112, 156, 114)))

# plot species colour-coded by respiration mode
(p <- ggtree(tree2, aes(color=group), layout = "fan", ladderize = TRUE, size=0.5) + 
    xlim(0,30) +
    scale_color_manual(values=c("cornflowerblue", "cornflowerblue", "chartreuse3"))+
    #scale_x_reverse() +
    geom_tiplab(size=2) +
    theme(legend.position = "none"))

# The tree object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

all_temp$phylo <- gsub(' ','_',all_temp$species_phylo)# On the tree, names have underscores instead of spaces.

## to compute branch length
all_temp.tree # inferred from taxonomic relatedness, no branch-lengths

is.binary(all_temp.tree) # check for polytomies

all_temp.tree_p <- multi2di(all_temp.tree) # resolve polytomies before calculating branch length

is.binary(all_temp.tree_p)

# plot new tree without polytomies
(all_temp.tree_p_2 <- groupClade(all_temp.tree_p, c(112, 156, 114)))

(p <- ggtree(all_temp.tree_p_2, aes(color=group), layout = "fan", ladderize = TRUE, size=0.5) + 
    xlim(0,30) +
    scale_color_manual(values=c("cornflowerblue", "cornflowerblue", "chartreuse3"))+
    #scale_x_reverse() +
    geom_tiplab(size=2) +
    theme(legend.position = "none"))

# compute branch length
all_temp.tree.bl <- compute.brlen(all_temp.tree_p, method="Grafen", power = 1) # Grafen method used to calculate branch-length from a non-ultrametric tree

# preliminary analyses showed that different tree transformation (rho) resulted in only slightly differences in model performance, so default is used (rho, 'power' = 1).

# This is how a transformation of rho = 1 looks like:
ggtree(all_temp.tree.bl, layout = "fan", open.angle=0, ladderize = TRUE, size=0.75) +
  geom_tiplab(size=2.1, hjust = -.1, color='blue')

# should we try a different tree transformation (rho different from 1) based on model performance, as in Verberk et al. 2022?

A <- vcv.phylo(all_temp.tree.bl) # calculate covariance matrix from tree

# explore whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.
setdiff(all_temp$phylo, rownames(A)) # check species names in dataset match those in A (0 is OK)
setdiff(rownames(A), all_temp$phylo) # check species names in A match those in dataset (0 is OK)

#-------------------------------------------------------------------------------

# Define priors

# mix of weakly informative and informative priors
priors =c(
  
  ## Informative priors:
  
  # Intercept:
  # Following MLBH (Glazier 2005, 2010), slope b approaches 2/3 with temperature-increased 
  # metabolic levels (log10 L = 0 (intercept), i.e. L of 1 mg O2/g/h, which is high):
  
  prior(normal(0.67, 0.5), "Intercept"),
  
  # Coefficient of log-transformed Metabolic level (L):
  # Following empirical estimate of slope b vs. L between species
  # at different temperatures from Killen et al. (2010) in Eco. Lett. (estimate = -0.145):
  
  prior(normal(-0.1, 0.5), "b", coef = "log10_L"),
  
  ## weakly informative priors
  prior(normal(0, 1), "b"),
  prior(student_t(3, 0, 2.5), "sd"),
  prior(student_t(3, 0, 2.5), "sigma")
)

#-------------------------------------------------------------------------------

#### Model fitting ####
if(fitt){
  # model 
  m_all.temp_complex.info <- brm(
    b ~ log10_L * group + (log10_L|experiment) + (1|gr(phylo, cov = A)),
    data = all_temp,
    family = student(),
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    data2 = list(A = A),
    prior = priors,
    sample_prior = TRUE,
    save_pars = save_pars(all = TRUE),
    iter = 7000, warmup = 3000, chains = 4, cores = 4)
  
  # total running time: 40 min # Bulk_ESS seems fine for all parameters (i.e. > 1000)
  
  # Model summary ####
  summary(m_all.temp_complex.info, prob = 0.95) %>% print(digits = 3)
  
  #  Family: student 
  #  Links: mu = identity; sigma = identity; nu = identity 
  #  Formula: b ~ log10_L * group + (log10_L | experiment) + (1 | gr(phylo, cov = A)) 
  #  Data: all_temp (Number of observations: 523) 
  #  Draws: 4 chains, each with iter = 7000; warmup = 3000; thin = 1;
  #  total post-warmup draws = 16000
  
  #  Group-Level Effects: 
  #    ~experiment (Number of levels: 149) 
  #                         Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #  sd(Intercept)             0.122     0.017    0.091    0.157 1.001     2590     5642
  #  sd(log10_L)               0.096     0.019    0.061    0.135 1.003     1073     2461
  #  cor(Intercept,log10_L)    0.415     0.193   -0.048    0.705 1.002     1302     2571
  
  #  ~phylo (Number of levels: 111) 
  #                Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #  sd(Intercept)    0.065     0.041    0.003    0.153 1.003     1293     3731
  
  #  Population-Level Effects: 
  #                              Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #  Intercept                      0.725     0.070    0.581    0.873 1.001     8649     6473
  #  log10_L                       -0.002     0.021   -0.043    0.038 1.001     6182     8591
  #  groupwaterMbreather           -0.027     0.085   -0.198    0.160 1.001     8600     6621
  #  log10_L:groupwaterMbreather   -0.090     0.028   -0.145   -0.035 1.001     5883     8162
  
  #  Family Specific Parameters: 
  #    Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #  sigma    0.031     0.004    0.024    0.039 1.001     2161     5709
  #  nu       1.782     0.266    1.339    2.374 1.000     4237     8438
  
  #  Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
  #  and Tail_ESS are effective sample size measures, and Rhat is the potential
  #  scale reduction factor on split chains (at convergence, Rhat = 1).
  
  # create a folder named 'Model_outputs', and save model (so that they don't need to be rerun for every session) ####
  # note: saved models are too large to be pushed in git standard repositories
  
  saveRDS(m_all.temp_complex.info,"../Model_outputs/m_all.temp.rds")
  m_all.temp <- readRDS("../Model_outputs/m_all.temp.rds") # rename saved model
  
  summary(m_all.temp, prob = 0.95) %>% print(digits = 3)
  
  # write stan code for this model:
  stancode(m_all.temp)
  
  # Posterior predictive checks ####
  plot(m_all.temp, N = 5, ask = FALSE) # check convergence of chains / looks fine (note sd_phylo_intercept)
  
  pp_check(m_all.temp, type="scatter_avg", ndraws = 100) # some outliers at the x end       
  pp_check(m_all.temp, type = "loo_pit_overlay", nsamples = NULL) # looks fine                             
  pp_check(m_all.temp, type = "dens_overlay", nsamples = 99) # looks fine
  
  plot(conditional_effects(m_all.temp), points = TRUE) # see data from the same experiment (Fig. 2)
  
  ## LOO 
  loo_model_all_temp <- loo(m_all.temp, save_psis = TRUE)
  loo_model_all_temp # 3 problematic observations (0.6% of total data) based on pareto-k-diagnostic (k > 0.7)
  
  plot(loo_model_all_temp) # plot of the pareto-k values
  
  ## RELOO without problematic observations, using one more chain:
  # (takes about 40 minutes to run, and results are almost identical)
  # This is a necessary step to calculate Monte Carlo SE of elpd_loo
  
  reloo_model_all_temp <- reloo(loo_model_all_temp, m_all.temp, chains = 1)
  reloo_model_all_temp # negligible improve of LOOIC (-1173.5 vs. -1173.8). 
  
  # Computed from 16000 by 523 log-likelihood matrix
  
  #           Estimate   SE
  #  elpd_loo    586.9 27.7
  #  p_loo       290.2 13.1
  #  looic     -1173.8 55.4
  #  ------
  #    Monte Carlo SE of elpd_loo is 0.3.
  
  #  Pareto k diagnostic values:
  #                           Count Pct.    Min. n_eff
  #  (-Inf, 0.5]   (good)     510   97.5%   197       
  #   (0.5, 0.7]   (ok)        13    2.5%   967       
  #     (0.7, 1]   (bad)        0    0.0%   <NA>      
  #     (1, Inf)   (very bad)   0    0.0%   <NA>      
  
  #  All Pareto k estimates are ok (k < 0.7).
  #  See help('pareto-k-diagnostic') for details.
  
}

if(fittsimulated){
  m_all.temp <- readRDS("../Model_outputs/m_all.temp.rds")
  trueparms <- as.list(posterior_summary(m_all.temp, variable = c("b_Intercept", "b_log10_L", "b_groupwaterMbreather", "b_log10_L:groupwaterMbreather", "sigma", "nu", "sd_experiment__Intercept", "sd_experiment__log10_L", "cor_experiment__Intercept__log10_L", "sd_phylo__Intercept"))[, 1])
  iter <- 10 #number of simulated data sets to generate
  simb <- array(dim = c(dim(all_temp)[1], iter))
  for(i in 1:iter){
    simdata <- simulatetdata(realdata = all_temp, trueparms = trueparms, A = A)
    simb[, i] <- simdata$b
    simfit <- brm(
      b ~ log10_L * group + (log10_L|experiment) + (1|gr(phylo, cov = A)),
      data = all_temp,
      family = student(),
      control = list(adapt_delta = 0.999, max_treedepth = 15),
      data2 = list(A = A),
      prior = priors,
      sample_prior = TRUE,
      save_pars = save_pars(all = TRUE),
      iter = 7000, warmup = 3000, chains = 4, cores = 4) # SAME NUMBER OF ITERATIONS THAN REAL MODEL
    saveRDS(simfit, paste(simfitname, "_", i, ".rds", sep = ""))
    
  }
  matplot(all_temp$b, simb, xlab = "real b", ylab = "simulated b", pch = 16) #THIS IS JUST FOR TESTING: ARE WE GENERATING PLAUSIBLE VALUES?
  abline(h = c(2/3, 1), lty = 2, lwd = 3) # typical boundaries of b under the MLBH
  abline(h = c(0.2, 1.6), lty = 2, lwd = 2, col = "red") # minimum and maximum value from real data
}

if(processtsimulated){
  fnames <- list.files(path = "../Model_outputs/", pattern = "simfitt_")
  m_all.temp <- readRDS("../Model_outputs/m_all.temp.rds")
  trueparms <- as.list(posterior_summary(m_all.temp, variable = c("b_Intercept", "b_log10_L", "b_groupwaterMbreather", "b_log10_L:groupwaterMbreather", "sigma", "nu", "sd_experiment__Intercept", "sd_experiment__log10_L", "cor_experiment__Intercept__log10_L", "sd_phylo__Intercept"))[, 1])
  plotsimulated(path = "../Model_outputs/", fnames = fnames, trueparms = trueparms)
  dev.off()
}

if(tdiagnostics){
  m_all.temp <- readRDS(file = "../Model_outputs/m_all.temp.rds")
  postnu <- as_draws_matrix(m_all.temp, variable = "nu")
  priornu <- rgamma(n = dim(postnu)[1], shape = 2, rate = 0.1) #sample from gamma distribution with parameters as in Stan code
  priornu <- priornu[priornu > 1] #truncate as in Stan code
  plot(density(postnu), xlab = expression(nu), ylab = "density", main = "", lwd = 2, xlim = range(priornu))
  lines(density(priornu), lwd = 2, lty = "dashed") #we can see that there is information on nu in the data: posterior is much more concentrated than prior
  legend("topright", bty = "n", lty = c("solid", "dashed"), lwd = c(2, 2), legend = c("posterior", "prior"))
  
  postsigma <- as_draws_matrix(m_all.temp, variable = "sigma")
  priorsigma <- abs(rt(n = dim(postsigma)[1], df = 3) * 2.5 + 0) #sample from half non-standardized t-distribution with parameters as in Stan code
  plot(density(postsigma), xlab = expression(sigma), ylab = "density", main = "", lwd = 2, xlim = range(priorsigma))
  lines(density(priorsigma), lwd = 2, lty = "dashed") #can see there is information on sigma in the data: posterior much more concentrated than prior
  legend("topright", bty = "n", lty = c("solid", "dashed"), lwd = c(2, 2), legend = c("posterior", "prior"))
  
  resid <- residuals(m_all.temp) #these are summarized: column 1 is estimated residual
  qqplot(x = qt(ppoints(dim(resid)[1]), df = mean(postnu)), y = resid[, 1], xlab = "theoretical", ylab = "observed")
  qqline(y = resid[, 1], distribution = function(p) qt(p, df = mean(postnu)), probs = c(0.25, 0.75))
  
  par(mfrow = c(4, 4))
  resid <- residuals(m_all.temp, summary = FALSE) #one row for each draw from posterior
  for(i in 1:16){
    qqplot(x = qt(ppoints(dim(resid)[2]), df = mean(postnu)), y = resid[i, ], xlab = "theoretical", ylab = "observed")
    qqline(y = resid[i, ], distribution = function(p) qt(p, df = mean(postnu)), probs = c(0.25, 0.75))
  }
  
}

if(fitgaussian){
  # Student-t distribution vs. Gaussian distribution ####
  # Check that model performance is better using a Student-t distribution than a Gaussian distribution:
  m_all.temp_gaussian <- update(m_all.temp, family = gaussian())
  
  # total running time: 25 min # Bulk_ESS are fine (i.e.> 1000)
  
  # save model (so that they don't need to be rerun for every session)
  # saveRDS(m_all.temp_gaussian,"../Model_outputs/m_all.temp_gaussian.rds")
  # m_all.temp_gaussian <- readRDS("../Model_outputs/m_all.temp_gaussian.rds") # call saved model
  
  # compare models through LOOIC
  loo(m_all.temp, m_all.temp_gaussian)
  
  # Model comparisons:
  #                             elpd_diff se_diff
  # m_all.temp_complex.info       0.0       0.0  
  # m_all.temp_complex_gaussian -94.4      18.5  
  
}
#-------------------------------------------------------------------------------
# Save data on the R session and packages versions for reproducibility shake ####
sink("../R_session/Bayesian_models_temp_R_session.txt") # Create folder named 'R_session' 
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################