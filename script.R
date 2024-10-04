# This R script runs with a renv environment available in a dedicated Github Repository. Below are the commands to recreate that environment, in case yo
#renv::init()
#renv::install("jsonlite")
#renv::install("MASS")
#renv::snapshot()

renv::restore()
# Set the model's parameters a priori
rho           <- 0.9
sigma_epsilon <- 0.1
beta0         <-
beta1         <-
alpha1        <-
alpha2        <-

# Generate error terms:
set.seed(20241004)
library(MASS)
ui1_ui2 <- mvrnorm(
           n = 100000, 
           mu = c(0, 0), 
           Sigma = matrix(c(sigma_epsilon^2, rho*sigma_epsilon, rho*sigma_epsilon, 1), 
                          ncol = 2)
                  )