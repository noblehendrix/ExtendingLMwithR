############################
#  Extending linear models with R
#  UDEC Dept. Oceanographico
#  Concepcion, Chile
#  8 - 14 January 2017
#  Prof: Noble Hendrix
#  noblehendrix@gmail.com
#--------------------------
#  Día 3 - Modelos lineales generalizados
#--------------------------

#########  Part I Binomial Modeling #########



## Some simulation of binomial data and a binomial GLM

#link function to relate mu.real to probability of success, p uses
# the logit function  logit(p) = mu.real, p = inv.logit(mu.real)
logit.fcn <- function(x) log(x/(1-x))
inv.logit.fcn <- function(x) exp(x)/(1 + exp(x))
#nota buena - both of these are also in the package boot()

#simulate a relationship between the covariate X and the probability of success p
#coefficient values
beta.true<- c(-1, 2)
#ordered values of the X covariate
X.vals <- seq(-2, 2, length = 30)
#the systematic or linear part of the model
mu.vals <- beta.true[1] + beta.true[2]*X.vals
#the link function to get the systematic part in the range of (0,1)
p.vals <- inv.logit.fcn(mu.vals)

#plot these relationships
par(mfrow=c(1,2))
plot(X.vals, mu.vals, main="Original")
plot(X.vals, p.vals, pch=16, main="Logit", ylim = c(0,1))

#Simulate the binomial trials 
set.seed(123)
#number of trials per X value (note that this does not need to be the same, can be different for each value of X)
#number of X values
nobs <- 50
N.trials <- 30
X.obs <- rnorm(nobs)
#the systematic or linear part of the model
mu.true <- beta.true[1] + beta.true[2]*X.obs
#the link function to get the systematic part in the range of (0,1)
p.true <- inv.logit.fcn(mu.true)

#Simulate the observed data
#Method 1: write a loop over the i = 1, ...,nobs 
Y.obs <- vector(length = nobs)
for(i in 1:nobs){
  Y.obs[i] <- rbinom(n=1,size=N.trials,prob=p.true[i])
}
# Method 2: use sapply() for looping over the i = 1,...,nobs 
Y.obs <- sapply(1:nobs, function(i) rbinom(n=1,size=N.trials,prob=p.true[i]))

#library(boot) - this has the function inv.logit() the function (logit),

par(mfrow = c(1,1))
p.obs <- Y.obs/N.trials
plot(X.obs, p.obs, pch=16, col=2, xlab = "X values", ylab = "Probability of Success")
curve(inv.logit.fcn(beta.true[1] + beta.true[2]*x), from=-3, to=13, add=TRUE, lwd=2)
successes <- Y.obs
failures <- N.trials-Y.obs

#Two ways to fit the binomial glm - 
## First method for specifiying a binomial GLM: use the vector of successes and the vector of failures together
(glm1 <- glm(cbind(successes, failures)~X.obs, family=binomial))
## Second method - use the successes/total and then put the number of trials in as weights
(glm2 <- glm(p.obs~X.obs, family=binomial, weights=rep(N.trials, length = nobs) ) )

summary(glm1)
summary(glm2)

#test of significance
drop1(glm2, test="Chi")

#test of significance
anova(glm2, test = "Chisq")

#evaluate the residuals
plot(glm2) 

#make some predictions
plot(X.obs, p.obs, pch=16, col=4, xlab = "X values", ylab = "Probability of Success")
#Make some predictions and plot them -  
newD <- data.frame(X.obs = seq(from = -3, to = 3, by = 0.1))
#Predictions on the scale of the response - probability of Tb in interval (0,1)
Pred <- predict(glm1,newdata = newD, type = "response")
lines(newD$X.obs,Pred)
#Predictions on the scale of the linear predictor (-Inf, Inf) then transform
# with inv.logit
#also we can specify se.fit = TRUE to obtain SE's that can be used to construct 95% CI 
Pred.LM <- predict(glm1, newdata = newD, type = "link", se.fit = TRUE)
p.pred.mean <- inv.logit.fcn(Pred.LM$fit)
p.pred.lo <- inv.logit.fcn(Pred.LM$fit - 1.96*Pred.LM$se.fit)
p.pred.hi <- inv.logit.fcn(Pred.LM$fit + 1.96*Pred.LM$se.fit)
#plot them
lines(newD$X.obs, p.pred.mean, col = 2)
lines(newD$X.obs, p.pred.lo, col = 2, lty = 2)
lines(newD$X.obs, p.pred.hi, col = 2, lty = 2)


####### 
## Fit an actual data set
#Tuberculosis, a respiratory disease, in boar (el jabalí)
Boar <- read.table("Boar.txt", header = TRUE)
#view the relationship between length and Tb prevalence
plot(x=Boar$LengthCT, y = Boar$Tb,xlab="Length", ylab="Tb", pch = 16)

##################
#Un reto pequeño I
#fit the binomial model: Tb is a function of Length


#test significance of Length on prevalence of Tb



########  Part II Poisson and Negative Binomial Modeling #########

#Simulate a Poisson distributed set of observations

X1.p <- rnorm(70)
X2.p <- rnorm(70)
X3.p <- X1.p*X2.p

#the systematic or linear component
log.lambda <- 3 + 2*X1.p - 1.5*X2.p + 0.8*X3.p

#the link
lambda.p <- exp(log.lambda)

#the distribution
Y.P <- rpois(70, lambda = lambda.p)

#construct the data frame 
data.p<- data.frame(Y = Y.P, X1 = X1.p, X2 = X2.p, X3 = X3.p)

#fit the model
#this is a forward selection where we add the terms sequentially to build the model 
glmp1<- glm(Y ~ X1, data = data.p, family = poisson)
anova(glmp1, test = "Chisq")

glmp2 <- update(glmp1, .~. + X2)

anova(glmp1, glmp2, test = "Chisq")

##################
#Un reto pequeño II
#How do we fit the interaction term?


#How could we perform a backwards selection - that is evaluate the full model to see if we should drop any terms?

#Poisson and negative binomial modeling of road kill data

#Roadkill data
#Data and analysis obtained from Ch 7 in Zuur et al. (2009)

RK <- read.table("Roadkills.txt", header = TRUE)
names(RK)

plot(RK$D.PARK,RK$TOT.N,xlab="Distance to park",
     ylab="Road kills")

#Start with a simple model that just focuses on distance from park
M1<-glm(TOT.N~D.PARK,family=poisson,data=RK)
summary(M1)

#plot model predictions
MyData=data.frame(D.PARK=seq(from=0,to=25000,by=1000))
G<-predict(M1,newdata=MyData,type="link",se=T)
F<-exp(G$fit)
FSEUP<-exp(G$fit+1.96*G$se.fit)
FSELOW<-exp(G$fit-1.96*G$se.fit)
#confidence plots
lines(MyData$D.PARK,F,lty=1)
lines(MyData$D.PARK,FSEUP,lty=2)
lines(MyData$D.PARK,FSELOW,lty=2)

#we want to look at the residual deviance relative to the residual #degrees of freedom to see if the assumptions of Poisson are met
summary(M1)
#They are not equal, so we have some more work to do...


#Now lets transform and use the rest of the possible covariates
RK$SQ.POLIC<-sqrt(RK$POLIC)
RK$SQ.WATRES<-sqrt(RK$WAT.RES)
RK$SQ.URBAN<-sqrt(RK$URBAN)
RK$SQ.OLIVE<-sqrt(RK$OLIVE)
RK$SQ.LPROAD<-sqrt(RK$L.P.ROAD)
RK$SQ.SHRUB<-sqrt(RK$SHRUB)
RK$SQ.DWATCOUR<-sqrt(RK$D.WAT.COUR)

#build a second model with all of the transformed covariates
M2<-glm(TOT.N~OPEN.L+MONT.S+SQ.POLIC+
         SQ.SHRUB+SQ.WATRES+L.WAT.C+SQ.LPROAD+
         SQ.DWATCOUR+D.PARK,family=poisson,data=RK)
summary(M2)

drop1(M2,test="Chi")

M3 <- glm(TOT.N ~ MONT.S + SQ.POLIC + D.PARK +
          SQ.SHRUB + SQ.WATRES + L.WAT.C + SQ.LPROAD +
          SQ.DWATCOUR, family = poisson, data = RK)
anova(M2, M3, test = "Chi")

#different way of writing the same model - 
M3a <- update(M2, .~. - OPEN.L )
anova(M3, M3a, test = "Chi")

drop1(M3a, test = "Chi")

summary(M3a)
#what is the residual deviance and degrees of freedom?  


#Models for over-dispersion:
#Zuur et al. (2009) use the quasipoisson distribution and an F test

M4<- glm(TOT.N ~ OPEN.L + MONT.S + SQ.POLIC+
         SQ.SHRUB + SQ.WATRES + L.WAT.C + SQ.LPROAD+
         SQ.DWATCOUR + D.PARK, family = quasipoisson, data = RK)
summary(M4)
#Note the estimate of the Disperson parameter
drop1(M4,test="F")

#
M5<- glm(TOT.N ~ MONT.S + SQ.SHRUB + L.WAT.C+ D.PARK, family = quasipoisson, data = RK)
summary(M5)
drop1(M5,test="F")

M6 <- update(M5, .~. - SQ.SHRUB)
summary(M6)
drop1(M6,test="F")


#Another option is to fit the models using glm.nb() in the library MASS

library(MASS)
M7<-glm.nb(TOT.N~OPEN.L+MONT.S+SQ.POLIC+
         SQ.SHRUB+SQ.WATRES+L.WAT.C+SQ.LPROAD+
         SQ.DWATCOUR+D.PARK,link="log",data=RK)
drop1(M7, test = "Chi")


#Select a final model - all terms with p-values < 0.01
M8<-glm.nb(TOT.N~OPEN.L+D.PARK,link="log",data=RK)
summary(M8)
drop1(M8,test="Chi")
par(mfrow=c(2,2))
plot(M8)

#Fit the same model as a Poisson - note these models are nested in that one has the dispersion parameter set to 1 (Poisson) and the other is estimating the dispersion parameter (Neg Binomial)
M9 <- glm(TOT.N ~ OPEN.L + D.PARK, family = poisson, data = RK)

#likelihood ratio test for nested models - 
llhNB = logLik(M8)
llhPoisson  =logLik(M9)
d <- 2 * (llhNB - llhPoisson)
pval <- pchisq(as.numeric(d), df=1, lower.tail=FALSE)/2


