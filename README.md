# Source code
R source files for *Vosseler, Alexander (2016):* [*Bayesian model selection for unit root testing with multiple structural breaks*](https://www.sciencedirect.com/science/article/abs/pii/S0167947314002485) published in *Computational Statistics & Data Analysis (CSDA)*.

*bayes.urtest.R*: Main wrapper file for the unit root testing procedure

*rjbreak_ADF5/6.R*: Core algorithm as described in the paper, Reversible jump MCMC, which is called in bayes.urtest.R among others

*Unemployment_rates_data_analysis_Rev_Jump_unitroot.R*: Shows how unit root test can be applied to real time series

