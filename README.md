# Intraspecific response of metabolic scaling to temperature and activity in water- and air-breathing ectothermic vertebrates
This repository contains codes and data needed to reproduce the analyses and figures of the manuscript:

García-Gómez G., et al. (to be completed), (2023). Combining theoretical approaches to understanding the intraspecific variation in metabolic scaling: responses to temperature and activity differ between water- and air-breathing ectothermic vertebrates. To be submitted to Ecology Letters.

## Cite the repository as:
(citation here)

## Scripts
* **Model_temperature_test.R** This script contains all the code necessary to reproduce the results from the Bayesian model analysing the relationship between the scaling slope and the *temperature-increased metabolic level* within ectothermic vertebrate species. 
* **Model_activity_test.R** This script contains all the code necessary to reproduce the results from the Bayesian model analysing the relationship between the scaling slope and the *activity-increased metabolic level* within ectothermic vertebrate species.

The latter scripts provide the option to run just the model with student-t (as in the present manuscript) or with Gaussian errors, as:
> fitt <- # set TRUE to fit the model with student-t errors

> fitgaussian <- # set TRUE to fit the model with Gaussian errors (pre-set as FALSE)

(see lines 49-50 in Model_temperature_test.R, and 41-42 in Model_activity_test.R)

Additionally, Model_temperature_test.R provides the option to (1) simulate 10 datasets under the temperature-effect model, with posterior mean parameter values and the same structure as the real data, and fitting the model to these simulated data sets; (2) plot estimates of the model fitted simulated data; and (3) check model residuals as:
> fittsimulated <- # set TRUE to fit the student-t model to data simulated under posterior mean parameters estimated from real data (pre-set as FALSE)

> processtsimulated <- # set TRUE to plot estimates from student-t model fitted to simulated data (pre-set as FALSE)

> tdiagnostics <- # set TRUE to get diagnostics for the model with student-t errors (pre-set as FALSE)

(see lines 51-53 in Model_temperature_test.R)

Last, the number of iterations in these models can be changed in line 180 of Model_temperature_test.R and line 169 of Model_activity_test.R.

## Data inputs
There are two datasets used in the analysis, which are included in the appendix of manuscript as **Table S1** and **Table S2**. These datasets contain data on intra-specific metabolic scaling of water- and air-breathing ectothermic vertebrates. Table S1 is used to test the effect of temperature, mediated by metabolic level, on the scaling slopes. Table S2 is used to test the effect of activity level, mediated by metabolic level, on the scaling slopes.

* **Table S1**, file: **table_S1.csv** contains data on metabolic scaling regressions performed at various temperatures from experiments in inactive animals. These experiments measured mostly resting or routine metabolism, which may include some spontaneous activity. The columns are:
  * *group* whether the species is water- or air-breather
  * *species* name of the species as reported in the original study
  * *species_phylo* name of the species as recorded in Open Tree of Life       (https://tree.opentreeoflife.org) 
  * *experiment* label for the experiment from which a set of scaling regressions was produced. Each experiment contains at least two regressions measured at two temperatures, whereas other conditions (e.g., metabolic state of animals) remain the same  
  * *temp* experimental temperature (in degrees Celsius) at which the scaling regression was measured. If the originial study reported temperature ranges instead of a single temperature for each regression, the mean temperatures of the ranges were used
  * *L* metabolic level (in mg O<sub>2</sub> g<sup>-1</sup> h<sup>-1</sup>) 
  * *log10_L* log<sub>10</sub>-transformed metabolic level (in mg O<sub>2</sub> g<sup>-1</sup> h<sup>-1</sup>) 
  * *mass_mid.g* body mass (in g) at the geometric midpoint of the mass range of the scaling regression
  * *min_mass.g* minimum body mass (in g) of the mass range of the scaling regression
  * *max_mass.g* maximum body mass (in g) of the mass range of the scaling regression
  * *L_state_study* metabolic state of animals as reported in the original study ("NA" if this state is not reported)
  * *study_id* identification number of the original study from which the data were compiled (see word file with the reference list for Table S1 and S2)

* **Table S2**, file: **table_S2.csv** contains data on metabolic scaling regressions performed at various temperatures from experiments in animals during different activity levels. These experiments measured from low activity levels as resting to routine metabolism (i.e., none or some spontaneous activity) to high activity levels as active and maximum metabolism (i.e., sustained activity or until exhaustion) by forcing locomotion. The columns are as in Table S1, plus three additional columns:
  * *experiment2* label for a set of scaling regressions measured at a single temperature in the same experiment, containing at least two regressions measured at different activity levels, whereas other conditions (e.g., metabolic state of animals) remain the same
  * *comp_L* label for the minimum and maximum activity levels that were measured in an experiment (e.g., "rest_max" means a experiment included measurements on resting and maximum metabolism).
  * *L_state* metabolic state of animals in each scaling regression as reported in the experiments 

## Notes
All data processing was carried out in the R software version 4.0.2. The R folder contains the scripts to reproduce the statistical analyses and figures presented in this manuscript. To improve clarity, figures were slightly edited by the addition of side annotations and shapes.

## R packages
The R packages used for each R script are enlisted in the corresponding R session files.

## Licence
This repository was provided by the authors under the (X replace with type of license) License.

## Further information
In case of further questions, please contact: Guillermo García-Gómez, email: guillegar.gz@gmail.com
