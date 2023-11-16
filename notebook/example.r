## testing for using package only
rm(list = ls())

library(devtools)
install_github("timktsang/serodynamics")
library(serodynamics)

###### input the data set

data("inputdata")
data("flu_activity")

###### run the MCMC to estimate parameter of the model
mcmc_result <- sero_dynamics(inputdata,flu_activity, 20000,10000,1)

##### obtain the model estimate from fitted MCMC result
extract_mcmc_result <- output_model_estimate(mcmc_result)
