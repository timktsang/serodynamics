<!-- README.md is generated from README.Rmd. Please edit that file -->

# serodynamics

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

`serodynamics` is a Bayesian modeling framework to infer influenza virus
infection status, antibody dynamics, and individual infection risks from
serological data, by accounting for measurement error. This could
identifiy influenza infections by relaxing 4-fold rise rule, and
quantifies the contributions of age and pre-epidemic
hemagglutination-inhibiting (HAI) titers to infection risk.

While the package currently provides a set of fundamental functions, we
are actively working on expanding its capabilities with more advanced
tools for a comprehensive understanding of your results.

## Installation

1.  Install \[R\]\[r-project\]

2.  Install the development version of serodynamics from
    [GitHub](https://github.com/timktsang/serodynamics):

<!-- -->

    devtools::install_github("timktsang/serodynamics")
    library(serodynamics)

## Example

This is a basic example of how to load some serological data and fitting
the model using the MCMC framework.

    library(serodynamics)

    ## Load in a data set and flu activity data
    data("inputdata")
    data("flu_activity")

    ###### run the MCMC to estimate parameter of the model
    ###### in actual analysis, number of iteration is 200000, burnin is 100000, and thinning is 10
    mcmc_result <- sero_dynamics(inputdata,flu_activity, 2000,1000,1)

    ##### obtain the model estimate from fitted MCMC result
    extract_mcmc_result <- output_model_estimate(mcmc_result)
