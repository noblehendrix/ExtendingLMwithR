############################
#  Extending linear models with R
#  UDEC Dept. Oceanographico
#  Concepcion, Chile
#  8 - 14 January 2017
#  Prof: Noble Hendrix & Cole Monnahan
#  noblehendrix@gmail.com
#--------------------------
#  Día 1 - Inferencía
#--------------------------


###### Example of normal

#Generate a random variable
#to reproduce results
set.seed(123)   

#A single sample 

rv.norm<- rnorm(1, mean = 20, sd = 2)

#Calculate probability density of mean = 15 given random sample
den.norm<- dnorm(x = rv.norm, mean = 15, sd = 2)

#Calculate probability density of mean = 15...25 given random sample - i.e. the likelihood
norm.range<- 15:25

den.norm.range<- dnorm(x = rv.norm, mean = norm.range, sd = 2)

#plot the probability density values
plot(x=norm.range, y=den.norm.range,  type = 'l', xlab = "values", ylab = "Density")

#Multiple samples 

#generate multiple samples from the normal distribution
rv.norm2<- rnorm(20, mean = 20, sd = 2)
#look at histogram of samples
hist(rv.norm2)

den.norm2<- dnorm(x = rv.norm2, mean = 15, sd = 2)

#calculate the product of the densities over all samples
prod(den.norm2)

#calculate the log density instead
log.den.norm2<- dnorm(x = rv.norm2, mean = 15, sd = 2, log = T)

#sum the log density
sum(log.den.norm2)

#Calculate probability density of mean = 15...25 given random samples using log density
#calculate length of the range
range.length<- length(norm.range)

#calculate log density = log likelihood (LL)
norm2.LL<- vector(length = range.length)
for(i in 1:11){
	norm2.LL[i]<- sum( dnorm(x = rv.norm2, mean = norm.range[i], sd = 2, log = T) )
	}

#plot the log density values
plot(norm.range, norm2.LL, xlab = "values", ylab = "Density")

#what value maximizes the log likelihood?


##############
# Estimation

#### Method 1  - develop likelihood and fit using optim

#Step 1 - develop function that returns negative log likelihood 
NLL_Norm1 = function(Par, Data){
  # Parameters
  mu = Par
  sigma = 2
  # Log-likelihood
  LL_i = dnorm( Data$Y, mean = mu, sd = sigma, log=TRUE )
  NLL = -1 * sum(LL_i)
  return( NLL )
}

Data1 = list( 'Y'= rv.norm2 ) #NOTE DATA MUST BE A LIST

#look at what NLL_Norm1 is doing - calculating the neg log like for different values of mu - note true value of mean is 20
NLL_Norm1(Par = 10, Data1)
NLL_Norm1(Par = 15, Data1)

vals<- seq(10, 30, by = 2)
range(vals)
NLL.vect<- vector(length = length(vals))

for(i in 1:length(vals)){
	NLL.vect[i]<- NLL_Norm1(Par = vals[i], Data1)
}

#Plot the NLL values over different hypotheses about value of p
plot(vals, NLL.vect, type = 'l', xlab = "p", ylab = "NLL")

#what is MLE = Min NLL?
vals[ which(NLL.vect == min(NLL.vect) ) ]

#Step 2 - let's use optim() to estimate mu for us
Start<- 15
Est = optim( par=Start, fn=NLL_Norm1, Data=Data1)
# Estimated parameters
print( Est$par ) # Estimated parameter


#what about estimating both the mean and the stdev of the sample? 


#Step 1 - develop function that returns negative log likelihood 
NLL_Norm2 = function(Par, Data){
  # Parameters
  mu = Par[1]
  sigma = Par[2]
  # Log-likelihood
  LL_i = dnorm( Data$Y, mean = mu, sd = sigma, log=TRUE )
  NLL = -1 * sum(LL_i)
  return( NLL )
}

#Step 2 - estimate parameters using optim
Start2<- c(15, 2.5)
Est2 = optim( par=Start2, fn=NLL_Norm2, Data=Data1, hessian=TRUE )
# Estimated parameters
print( Est2$par ) # Estimated parameter
print( Est2$hessian ) #Estimated Hessian matrix

# Estimated standard errors
print( sqrt(diag( solve(Est2$hessian) )) ) # square root of diagonal elements of the inverse-hessian matrix


###################
#### Example of sizes from two populations

#simulate sizes
mean1 <- 45
mean2 <- 55
stdev <- 4

#plots of the underlying distributions
curve(dnorm(x, mean1, stdev), from = 30, to = 80,  col = 2)
curve(dnorm(x, mean2, stdev), from = 30, to = 80, add = TRUE, col = 4)

# samples from the distributions
nobs <- 100
pop1 <- rnorm(nobs, mean1, sd = stdev)
pop2 <- rnorm(nobs, mean2, sd = stdev)

hist(pop1, border = 2, breaks = seq(30,80, by = 2), prob = TRUE, ylim = c(0, 0.15), main = "Size Samples", xlab = "Size")
rug(pop1, col = 6)
hist(pop2, border = 4, breaks = seq(31,81, by = 2), prob = TRUE, add = TRUE)
rug(pop2, col = 7)
legend("topright", c("North", "South"), col = c(4,2), lwd = c(1,1))

#Construct a dataset 
Dat1 <- list(Y = pop1)
Dat2 <- list(Y = pop2)

#Estimate parameters of the model - mu1, mu2, sigma

#Step 1 - develop function that returns negative log likelihood 
NLL_Norm3 = function(Par, Data1, Data2){
  # Parameters
  mu1 = Par[1]
  mu2 = Par[2]
  sigma = Par[3]
  # Log-likelihoods
  LL1 = dnorm( Data1$Y, mean = mu1, sd = sigma, log=TRUE )
  LL2 = dnorm( Data2$Y, mean = mu2, sd = sigma, log=TRUE )
  
  NLL = -1 * sum(LL1, LL2)
  return( NLL )
}

#Step 2 - estimate parameters using optim
Start3<- c(40, 50, 3)
#check starting conditions
NLL_Norm3(Start3, Data1 = Dat1, Data2 = Dat2)

Est3 = optim( par=Start3, fn=NLL_Norm3, Data1=Dat1, Data2 = Dat2, hessian=TRUE)
# Estimated parameters
print( Est3$par ) # Estimated parameter
print( Est3$hessian ) #Estimated Hessian matrix
print( sqrt(diag( solve(Est3$hessian) )) ) # square root of diagonal elements of the inverse-hessian matrix

#########
# Compare two models 

#Fit the 2 population data using a single mean and st dev
all.data <- list(Y = c(pop1, pop2))
Start2<- c(40, 2.5)
Est3a = optim( par=Start2, fn=NLL_Norm2, Data=all.data, hessian=TRUE )
# Estimated parameters
print( Est3a$par ) # Estimated parameter
print( Est3a$hessian ) #Estimated Hessian matrix
print( sqrt(diag( solve(Est3a$hessian) )) )

#compare neg log likelihood values - 
Est3$val
Est3a$val

#compute the likelihood ratio statistic for the 2 models
# D = 2*(log likelihood complex - log like simple)
# we computed the negative log likelihood values in our functions
D.val <- 2*(-Est3$val - (-Est3a$val))
#compare to Chi-sq with df_alt - df_simple = 3 - 2 = 1 df
p.val <- pchisq(D.val, df = 1, lower.tail = FALSE)


#plot the sample and the distributions implied by the results
hist(pop1, border = 2, breaks = seq(30,80, by = 2), prob = TRUE, ylim = c(0, 0.15), main = '', xlab = "Size")
rug(pop1, col = 2)
hist(pop2, border = 4, breaks = seq(31,81, by = 2), prob = TRUE, add = TRUE)
rug(pop2, col = 4)
#legend("topright", c("North", "South"), col = c(4,2), lwd = c(1,1))

#plots of the underlying distributions
#different means
curve(dnorm(x, Est3$par[1], Est3$par[3]), from = 30, to = 80,  col = 2, add = TRUE, lwd = 2)
curve(dnorm(x, Est3$par[2], Est3$par[3]), from = 30, to = 80, add = TRUE, col = 4, lwd = 2)
#same mean
curve(dnorm(x, Est3a$par[1], Est3a$par[2]), from = 30, to = 80, add = TRUE, col = 1, lwd = 2)
legend("topright", c("North", "South", "Combined"), col = c(4,2,1), lwd = c(1,1,1))

