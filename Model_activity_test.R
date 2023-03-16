################################################################################
# Script to fit phylogenetic multilevel model using brms: Part 2
# Title: Activity effect on metabolic scaling of ectothermic vertebrates
# Authors: Guillermo García-Gómez (guillegar.gz@gmail.com) & Matthew Spencer
# Date: 16032023
# Operating System: Windows 10 Pro 21H2
# ------------------------------------------------------------------------------
# Cite as: ?
# ------------------------------------------------------------------------------
rm(list=ls())# clear the work environment
# ------------------------------------------------------------------------------
# It needs to be set to Project directory
getwd()# to check
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

fitt <- TRUE #fit model with student-t errors?
fitgaussian <- FALSE #fit model with Gaussian errors?

# ------------------------------------------------------------------------------
set.seed(1234)# we need this to replicate the results

# General STAN specifications
rstan::rstan_options(auto_write = TRUE) # write objects to the hard disk to avoid recompilation 
options(mc.cores = parallel::detectCores())

#### Load data ####
all_act <- read.csv("table_S2.csv") # Table S2

# Histograms: 

# response variable: slopes b
hist(all_act$b)

# predictor variable: log10(Metabolic level, L)
hist(log10(all_act$L))

#### Models for all species (water + air-breathers) ####

# Building the phylogenetic tree for all species (needs internet connection)

all_act.taxa <- tnrs_match_names(as.character(all_act$species_phylo)) # match species from dataset to those included in Tree of Life
all_act.tree <- tol_induced_subtree(ott_ids = ott_id(all_act.taxa), label = "name") # building a tree

# Check phylogenetic tree
ggtree(all_act.tree, layout = "fan", open.angle=0, ladderize = TRUE, size=0.75) +
  geom_tiplab(size=2.1, hjust = -.1, color='blue')

# group species colour-coded according to respiration mode by selecting nodes in the tree 
(tree_plot <- groupClade(all_act.tree, c(59, 50)))

# plot species colour-coded by respiration mode
(p <- ggtree(tree_plot, aes(color=group), layout = "fan", ladderize = TRUE, size=0.5) + 
    xlim(0,30) +
    scale_color_manual(values=c("cornflowerblue", "cornflowerblue", "chartreuse3"))+
    #scale_x_reverse() +
    geom_tiplab(size=2) +
    theme(legend.position = "none"))

# The tree object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

all_act$phylo <- gsub(' ','_',all_act$species_phylo)# On the tree, names have underscores instead of spaces.

## to compute branch length
all_act.tree # inferred from taxonomic relatedness, no branch-lengths

is.binary(all_act.tree) # check for polytomies

all_act.tree_p <- multi2di(all_act.tree) # resolve polytomies before calculating branch length

is.binary(all_act.tree_p)

# plot new tree without polytomies
(all_act.tree_plot2 <- groupClade(all_act.tree, c(59, 50)))

(p2 <- ggtree(all_act.tree_plot2, aes(color=group), layout = "fan", ladderize = TRUE, size=0.5) + 
    xlim(0,30) +
    scale_color_manual(values=c("cornflowerblue", "cornflowerblue", "chartreuse3"))+
    #scale_x_reverse() +
    geom_tiplab(size=2) +
    theme(legend.position = "none"))

# compute branch length
all_act.tree.bl <- compute.brlen(all_act.tree_p, method="Grafen", power = 1) # Grafen method used to calculate branch-length from a non-ultrametric tree

# preliminary analyses showed that different tree transformation (rho) resulted in only slightly differences in model performance, so default is used (rho, 'power' = 1).

# This is how a transformation of rho = 1 looks like:
ggtree(all_act.tree.bl, layout = "fan", open.angle=0, ladderize = TRUE, size=0.75) +
  geom_tiplab(size=2.1, hjust = -.1, color='blue')

# should we try a different tree transformation (rho different from 1) based on model performance, as in Verberk et al. 2022?

B <- vcv.phylo(all_act.tree.bl) # covariance matrix from tree

# explore whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.
setdiff(all_act$phylo, rownames(B)) # check species names in dataset match those in A (0 is OK)
setdiff(rownames(B), all_act$phylo) # check species names in A match those in dataset (0 is OK)

#-------------------------------------------------------------------------------

# Define priors

# mix of weakly informative and informative priors
priors =c(
  
  ## Informative priors:
  
  # Intercept:
  # Following MLBH (Glazier 2005, 2010), slope b approaches 2/3 with temperature-increased 
  # metabolic levels (log10 L = 0 (intercept), i.e. L of 1 mg O2/g/h, which is high):
  
  prior(normal(1, 0.5), "Intercept"),
  
  # Coefficient of log-transformed Metabolic level (L):
  # Following empirical estimate of slope b vs. L between species
  # at different temperatures from Killen et al. (2010) in Eco. Lett. (estimate = -0.145),
  # we apply the inverse of this value for activity:
  
  prior(normal(0.1, 0.5), "b", coef = "log10_L"),
  
  ## weakly informative priors
  prior(normal(0, 1), "b"),
  prior(student_t(3, 0, 2.5), "sd"),
  prior(student_t(3, 0, 2.5), "sigma")
)

#-------------------------------------------------------------------------------

#### Model fitting ####
if(fitt){
  # model 
  m_all.act_complex.info <- brm(
    b ~ log10_L * group + temp + (log10_L|experiment) + (1|gr(phylo, cov = B)),
    data = all_act,
    family = student(),
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    data2 = list(B = B),
    prior = priors,
    sample_prior = TRUE,
    save_pars = save_pars(all = TRUE),
    iter = 7000, warmup = 3000, chains = 4, cores = 4)
  
  # total running time: 15 min # Bulk_ESS seems fine for all parameters (i.e. > 1000) 
  
  # Model summary ####
  summary(m_all.act_complex.info, prob = 0.95) %>% print(digits = 3)
  
  #  Family: student 
  #  Links: mu = identity; sigma = identity; nu = identity 
  #  Formula: b ~ log10_L * group + temp + (log10_L | experiment) + (1 | gr(phylo, cov = B)) 
  #  Data: all_act (Number of observations: 281) 
  #  Draws: 4 chains, each with iter = 7000; warmup = 3000; thin = 1;
  #  total post-warmup draws = 16000
  
  #  Group-Level Effects: 
  #    ~experiment (Number of levels: 56) 
  #                         Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #  sd(Intercept)             0.103     0.028    0.044    0.157 1.003     1804     2244
  #  sd(log10_L)               0.141     0.026    0.092    0.194 1.001     3438     4153
  #  cor(Intercept,log10_L)    0.827     0.127    0.567    0.966 1.001     2387     3009
  
  #  ~phylo (Number of levels: 47) 
  #                 Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #   sd(Intercept)    0.136     0.053    0.041    0.247 1.003     1367     2949
  
  #  Population-Level Effects: 
  #                              Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #  Intercept                      0.978     0.145    0.679    1.273 1.000    10953     8947
  #  log10_L                        0.194     0.055    0.085    0.302 1.001     6758     8638
  #  groupwaterMbreather           -0.111     0.162   -0.438    0.225 1.000     8481     7272
  #  temp                          -0.002     0.001   -0.004    0.001 1.001    11530    12005
  #  log10_L:groupwaterMbreather   -0.152     0.062   -0.274   -0.029 1.000     6639     9717
  
  #  Family Specific Parameters: 
  #        Estimate Est.Error l-95% CI u-95% CI  Rhat Bulk_ESS Tail_ESS
  #  sigma    0.059     0.007    0.046    0.073 1.001     5281     8858
  #  nu       4.777     2.318    2.358   10.667 1.000     9518    10820
  
  #  Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
  #  and Tail_ESS are effective sample size measures, and Rhat is the potential
  #  scale reduction factor on split chains (at convergence, Rhat = 1).
  
  # create a folder named 'Model_outputs', and save model (so that they don't need to be rerun for every session) ####
  # note: saved models are too large to be pushed in git standard repositories, save them elsewhere
  
  saveRDS(m_all.act_complex.info,"../Model_outputs/m_all.act.rds")
  m_all.act <- readRDS("../Model_outputs/m_all.act.rds") # rename saved model
  
  summary(m_all.act, prob = 0.95) %>% print(digits = 3)
  
  # write stan code for this model:
  stancode(m_all.act)
  
  # Posterior predictive checks ####
  plot(m_all.act, N = 6, ask = FALSE) # check convergence of chains / looks fine (but note sd_phylo_intercept)
  
  pp_check(m_all.act, type="scatter_avg", ndraws = 100) # some outliers at the x end       
  pp_check(m_all.act, type = "loo_pit_overlay", nsamples = NULL) # looks fine                             
  pp_check(m_all.act, type = "dens_overlay", nsamples = 99) # looks fine
  
  plot(conditional_effects(m_all.act), points = TRUE) # see data from the same experiments (Fig. 2)
  
  ## LOO 
  loo_model_all_act <- loo(m_all.act, save_psis = TRUE)
  loo_model_all_act # 3 problematic observations (1.1% of total data) based on pareto-k-diagnostic (k > 0.7)
  
  plot(loo_model_all_act) # plot of the pareto-k values
  
  ## RELOO without problematic observations, using one more chain:
  # (takes about 20 minutes to run, and results are almost identical)
  # This is a necessary step to calculate Monte Carlo SE of elpd_loo
  
  reloo_model_all_act <- reloo(loo_model_all_act, m_all.act, chains = 1)
  reloo_model_all_act # almost identical LOOIC (-560.4 vs. -560.6). 
  
  #  Computed from 16000 by 281 log-likelihood matrix
  
  #           Estimate   SE
  #  elpd_loo    280.3 16.3
  #  p_loo        77.3  5.6
  #  looic      -560.6 32.6
  #  ------
  #    Monte Carlo SE of elpd_loo is 0.2.
  
  #  Pareto k diagnostic values:
  #                           Count Pct.    Min. n_eff
  #    (-Inf, 0.5]   (good)     267   95.0%   242       
  #     (0.5, 0.7]   (ok)        14    5.0%   252       
  #       (0.7, 1]   (bad)        0    0.0%   <NA>      
  #       (1, Inf)   (very bad)   0    0.0%   <NA>      
  
  # All Pareto k estimates are ok (k < 0.7).
  # See help('pareto-k-diagnostic') for details.
}

if(fitgaussian){
  # Student-t distribution vs. Gaussian distribution ####
  # Check that model performance is better using a Student-t distribution than a Gaussian distribution:
  m_all.act_gaussian <- update(m_all.act, family = gaussian())
  
  # total running time: 8 min # Bulk_ESS are fine (i.e.> 1000)
  
  # save model (so that they don't need to be rerun for every session)
  # saveRDS(m_all.act_gaussian,"../Model_outputs/m_all.act_gaussian.rds")
  # m_all.act_gaussian <- readRDS("../Model_outputs/m_all.act_gaussian.rds") # call saved model
  
  # compare models through LOOIC
  loo(m_all.act, m_all.act_gaussian)
  
  # Model comparisons:
  #                   elpd_diff se_diff
  # m_all.act            0.0       0.0  
  # m_all.act_gaussian -15.9       7.0  
  
}


#-------------------------------------------------------------------------------
# Save data on the R session and packages versions for reproducibility shake ####
sink("../R_session/Bayesian_models_act_R_session.txt") # Create folder named 'R_session' 
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################