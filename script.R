# This R script runs with a renv environment available in a dedicated Github Repository. Below are the commands to recreate that environment, in case yo
#install.packages(renv)
#renv::init()
#renv::install("languageserver")
#renv::install("httpgd")
#renv::install("jsonlite")
#renv::install("MASS")
#renv::snapshot()


renv::restore()

################################################
# (a)

# Set the model's parameters a priori
rho           <- 0.9
sigma_epsilon <- 0.1
beta0         <- 4
beta1         <- 5
alpha0        <- -12 # Will be tuned later.
alpha1        <- 1
alpha2        <- 1

# Generate Ei
set.seed(202410041)
Ei <- round(abs(rnorm(100000, 9, 5)))

# Generate Zi
set.seed(202410042)
Zi <- rnorm(100000, 3, 0.75)

# Generate error terms:
set.seed(202410043)
library(MASS)
ui1_ui2 <- mvrnorm(
           n = 100000, 
           mu = c(0, 0), 
           Sigma = matrix(c(sigma_epsilon^2, rho*sigma_epsilon, rho*sigma_epsilon, 1), 
                          ncol = 2)
                  )

# Compute Wi
Wi <- beta0+beta1*Ei+ui1_ui2[,1]

# Fine tune alpha0 such that the participation condition is met by exactly half of the data
participation_condition <- alpha0+alpha1*Ei+alpha2*Zi+ui1_ui2[,2]
length(participation_condition[participation_condition>0])
# We can see that, with this value of alpha0 there are too many observations that satisfy the participation condition.
while (length(participation_condition[participation_condition>0])>50000) {
   alpha0 <- alpha0-0.00001
   participation_condition <- alpha0+alpha1*Ei+alpha2*Zi+ui1_ui2[,2]
   print(alpha0)
   print(length(participation_condition[participation_condition>0]))
} 
# We get alpha0=-12.03036. Explicitely set alpha0 accordingly:
alpha0 <- -12.03037
participation_condition <- alpha0+alpha1*Ei+alpha2*Zi+ui1_ui2[,2]
length(participation_condition[participation_condition>0])


# Create Wi_observed and participate dummy according to the participation condition
Wi_observed <- Wi
participate <- Wi/Wi
for (i in 1:length(participation_condition))
{
  if (participation_condition[i]<=0) {Wi_observed[i] <-NA}
  if (participation_condition[i]<=0) {participate[i] <-0}
}

# Computes the inverse Mills ratio for the participation equation
standard_normal_probability_density <- function(x) {
    (1/sqrt(2*pi))*exp(-(x^2)/(2))
}
inverse_mills_ratio <- dnorm((alpha0+alpha1*Ei+alpha2*Zi))/(pnorm((alpha0+alpha1*Ei+alpha2*Zi))

# Create the dataset
dataset <- data.frame(Wi=Wi,Wi_observed=Wi_observed,participate=participate,Ei=Ei,Zi=Zi,ui1=ui1_ui2[,1],ui2=ui1_ui2[,2],participation_condition=participation_condition,inverse_mills_ratio=inverse_mills_ratio)

# Estimate the parameters of the three regressions
true_model <- lm(Wi~Ei, data=dataset)
summary(true_model)
model_those_participating <- lm(Wi_observed~Ei, data=dataset)
summary(model_those_participating)
model_those_participating_control <- lm(Wi_observed~Ei+inverse_mills_ratio, data=dataset)
summary(model_those_participating_control)

# Interpret the regression parameters and the direction of the bias

################################################
# (b)

# Estimate with probit the participation decision model:
participation_decision <- glm(participate~Ei+Zi,data = dataset,family = binomial(link = probit))
participation_probability <- predict(participation_decision,dataset,type="response")
dataset <- data.frame(Wi=Wi,Wi_observed=Wi_observed,participate=participate,Ei=Ei,Zi=Zi,ui1=ui1_ui2[,1],ui2=ui1_ui2[,2],participation_condition=participation_condition,inverse_mills_ratio=inverse_mills_ratio, participation_probability=participation_probability)

# We can draw the plots:
plot(inverse_mills_ratio~participation_probability, data=dataset)
plot(inverse_mills_ratio~participation_probability, data=subset(dataset, alpha0+alpha1*Ei+alpha2*Zi+ui1_ui2[,2]>0))

# Discuss the non-linearity of the inverse Mills ratio.

################################################
# (c)

# The best fit line for the sample who participate:
best_fit_line_participants <- lm(inverse_mills_ratio~participation_probability, data=subset(dataset, alpha0+alpha1*Ei+alpha2*Zi+ui1_ui2[,2]>0))
plot(inverse_mills_ratio~participation_probability, data=subset(dataset, alpha0+alpha1*Ei+alpha2*Zi+ui1_ui2[,2]>0))
plot(dataset[inverse_mills_ratio],dataset[participation_probability])
abline(a=coef(best_fit_line_participants)[1], b=coef(best_fit_line_participants)[2])

abline(best_fit_line_participants)
plot(best_fit_line_participants)
plot(inverse_mills_ratio~participation_probability, data=dataset)
abline(dataset[inverse_mills_ratio]~dataset[participation_probability])

library(ggplot2)

ggplot(dataset, aes(x=participation_probability, y=inverse_mills_ratio)) +
    geom_point() + geom_smooth(method=lm, se=FALSE)

# The best fit line for the full sample.
best_fit_line_full_sample <- lm(inverse_mills_ratio~participation_probability, data=dataset)
plot(inverse_mills_ratio~participation_probability, data=dataset)
abline(a=coef(best_fit_line_full_sample)[1], b=coef(best_fit_line_full_sample)[2])


#Write a brief description.

################################################
# (d)

# Attach all of your code.