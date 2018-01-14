############################
#  Extending linear models with R
#  UDEC Dept. Oceanographico
#  Concepcion, Chile
#  8 - 14 January 2017
#  Prof: Noble Hendrix
#  noblehendrix@gmail.com
#--------------------------
#  Día 4 - Effectos aleatorios
#--------------------------

#simulate data and estimate parameters
#we will use lme4 for this example 
library(lme4)
set.seed(123)
#Generate some lengths - 30 for each population
length <- runif(30, 0, 20) 
#Generate 3 random intercepts with a mean of 40 and sd of 10
b0<-rnorm(3, mean = 40, sd = 10)
#Use the same slope for all populations
b1<- 2.0
#Calculate the true mass for each population
mass.true <- sapply(1:3, function(x){b0[x]+b1*length})
mass.obs <- mass.true + matrix(nrow = 30, ncol = 3, data = rnorm(90, mean = 0, sd = 10) )

plot(length, mass.obs[,1], col = 2, xlim = c(-1, 20), ylim = c(0, 120), cex = 1.5, lwd = 2, pch = 16, xlab = "Length", ylab = "Mass")
points(length, mass.obs[,2], col = 3, pch = 16,cex = 1.5)
points(length, mass.obs[,3], col = 4, pch = 16,cex = 1.5)


#Make the matrix of mass values into a vector
mass.obs.vec<- c(mass.obs)
length.vec <- rep(length, times = 3)
pop <- factor(rep(1:3, each = 30) )

sim.data <- data.frame(mass = mass.obs.vec, length = length.vec, pop = pop )

# Fit mixed model, print random effects and plot regression lines
lmm <- lmer(mass ~ length + (1|pop), data = sim.data)
#compare estimates of intercept to true values
b0  #true values
fixef(lmm)[1] + ranef(lmm)$pop   #estimates

abline((fixef(lmm)[1]+ranef(lmm)$pop)[1,], fixef(lmm)[2], col = 2, lwd = 3)
abline((fixef(lmm)[1]+ranef(lmm)$pop)[2,], fixef(lmm)[2], col = 3, lwd = 3)
abline((fixef(lmm)[1]+ranef(lmm)$pop)[3,], fixef(lmm)[2], col = 4, lwd = 3)

#overall mean
abline(lm(mass~length, data = sim.data) , lwd = 3, lty = 3)

#other designs for RE
#random slope on length only
#(0+length|pop)
#both intercept and length
#(1+length|pop)

#The following example was from:
#    Mixed Effects Models and Extensions in Ecology with R (2009)
#    Zuur, Ieno, Walker, Saveliev and Smith.    Springer
#    This file was produced by Alain Zuur (highstat@highstat.com)
#    www.highstat.com


#we will use nlme for this example 
RIKZ <- read.table("RIKZ.txt", header = TRUE)

#Exposure is a per site variable, 

library(nlme)
RIKZ$fBeach <- factor(RIKZ$Beach)
RIKZ$fExp<-RIKZ$Exposure
RIKZ$fExp[RIKZ$fExp==8]<-10
RIKZ$fExp<-factor(RIKZ$fExp,levels=c(10,11))

# Richness - species richness (dependent variable)
# NAP - the height of a sampling station compared to mean tidal level
# Exposure - amount of exposure per beach
# Beach - beach number
# fExp - exposure as a factor of 2 levels
# 
#Step 1 - fit a full model with all interactions (note the syntax 
#of NAP*fExp - NAP + fExp + NAP:fExp) of covariates and no random effects 
#using REML- use gls in the nlme package since no random effects 
B1 <- gls(Richness ~ 1 + NAP * fExp,
          method = "REML", data = RIKZ)
#Step 2a - fit a model with all interactions of covariates and the first random effect structure with REML
B2 <- lme(Richness ~1 + NAP * fExp, data = RIKZ,
          random = ~1 | fBeach, method = "REML")
#Step 2b - fit a model with all interactions of covariates and the second random effect structure with REML
B3 <- lme(Richness ~ 1 + NAP * fExp,data = RIKZ,
          random = ~1 + NAP | fBeach, method = "REML")

#Step 3a - evaluate the different models with REML using AIC and anova
AIC(B1,B2,B3)
anova(B1,B2,B3)

#Step 3b -identify the random effect structure of model B2

#Step 4 - switch to using ML and evaluate the covariate structure in the linear part of the model
#modify some parameters of the lme fitting to increase iterations:
lmc <- lmeControl(niterEM = 5200, msMaxIter = 5200)
#refit model B2 above using method = ML
M4A <- lme(Richness ~1 + NAP * fExp, control = lmc, data = RIKZ,
           random = ~1 | fBeach, method = "ML")

#use summary and anova to evaluate the models
summary(M4A)
anova(M4A)

#perhaps remove the interaction effect
M4B <- lme(Richness ~1 + NAP + fExp, control = lmc, data = RIKZ,
           random = ~1 | fBeach, method = "ML")

anova(M4A, M4B)

#Step 5 - identify the final model structure with respect to the covariates, Zuur selects model M4B

#Step 6 - rerun the final model with method = REML

M6 <- lme(Richness ~1 + NAP + fExp, control = lmc, data = RIKZ,
                 random = ~1 | fBeach, method = "REML")
#a plot of the standardized residuals
plot(M6)


###########
####
# GLMM - generalized linear mixed model

data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014, sd = 0.1){
  # n: Number of years
  # alpha, beta1, beta2, beta3: coefficients of a 
  #    cubic polynomial of count on year
  # sd: standard deviation of normal distribution assumed for year effects
  
  # Generate values of time covariate
  year <- 1:n
  
  # First level of noise: generate random year effects
  eps <- rnorm(n = n, mean = 0, sd = sd)
  
  # Signal (plus first level of noise): build up systematic part of the GLM and add the random year effects
  log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3 + eps
  expected.count <- exp(log.expected.count)
  mean.count <- exp(alpha + beta1 * year + beta2 * year^2 + beta3 * year^3 )
  # Second level of noise: generate random part of the GLM: Poisson noise around expected counts
  C <- rpois(n = n, lambda = expected.count)
  
  # Plot simulated data
  plot(year, C, type = "b", lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year", ylim = c(0, 1.1*max(C)))
  lines(year, expected.count, type = "l", lwd = 3, col = "red")
  lines(year, mean.count, col = "green", lwd=3)
  return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, sd = sd, expected.count = expected.count, C = C))
}

data <- data.fn()

library(lme4)
yr <- factor(data$year)         # Create a factor year
glmm.fit <- glmer(C ~ (1 | yr) + year + I(year^2) + I(year^3), family = poisson, data = data)

#success?

mny <- mean(data$year)
sdy <- sd(data$year)
cov1 <- (data$year - mny) / sdy
cov2 <- cov1 * cov1
cov3 <- cov1 * cov1 * cov1
glmm.fit <- glmer(C ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson, data = data)
glmm.fit

# Plot simulated data
plot(data$year, data$C, type = "b", lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year", ylim = c(0, 1.1*max(data$C)))
lines(data$year, data$expected.count, type = "l", lwd = 3, col = "red")

R.predictions <- exp(fixef(glmm.fit)[1] + fixef(glmm.fit)[2]*cov1 + fixef(glmm.fit)[3]*cov2 + fixef(glmm.fit)[4]*cov3 + unlist(ranef(glmm.fit)))
lines(data$year, R.predictions, col = "green", lwd = 2, type = "l")


#Fit to the actual falcons data set - 

peregrine <- read.table("Falcons.txt", header = TRUE)
yr <- factor(peregrine$Year)
mny <- mean(peregrine$Year)
sdy <- sd(peregrine$Year)
#construct a scaled covariate: subtract the mean and divide by the sd
cov1 <- (peregrine$Year - mny) / sdy
cov2 <- cov1 * cov1
cov3 <- cov1 * cov1 * cov1


################
## El Reto o Desafio
#What do we do next to model the actual data?
data.hal<-data.frame(Y = peregrine$Pairs, X1 = cov1, X2 = cov2, X3 = cov3)
#How do we evaluate the covariates?
# glmm.fit <- glmer(C ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson, data = data)
glmm.hal<- glmer(Y~  (1|yr) + X1 + X2 + X3, family = poisson , data = data.hal)

#modelo sin X3
glmm.hal1<- update(glmm.hal, .~. - X3)
anova(glmm.hal, glmm.hal1)
