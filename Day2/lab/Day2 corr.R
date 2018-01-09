############################
#  Extending linear models with R
#  UDEC Dept. Oceanographico
#  Concepcion, Chile
#  8 - 14 January 2017
#  Prof: Noble Hendrix
#  noblehendrix@gmail.com
#--------------------------
#  Día 2 - Modelos lineales
#--------------------------


#-------   Part I   ----------
## Simple linear regression

#simulate data
x <- runif(50, 0, 10)
y <- 1.5+3*x+rnorm(50)

#plot the simulated data
plot(x,y)

#a way to explore the linear model - thanks to Cole:
plot.example <- function(slope, intercept){
  plot(x,y)
  errors <- (intercept + x*slope) - y
  SSRes <- sum(errors^2)
  abline(a=intercept, b=slope)
  text(x=-20, y=20, labels=paste("Sum Sq Resids=", round(SSRes,2)))
}

plot.example(intercept = 3, slope = 1.5)


library(manipulate)
manipulate(plot.example(slope, intercept), slope=slider(0,10),
           intercept=slider(-5,5))

#fit a linear model
lm1<- lm(y~x)

#You can add a linear fitted Line to the existing scatter plot
abline(lm1, col = 2)
lm1
#force intercept through the origin
lm2<- lm(y ~ -1 + x)

#add this to the plot too
abline(lm2, lty = 2, col = 3)

## Add a legend to identify the two different lines
legend("bottomright",col=c(2,3), lty=1:2,
       legend=c("linear model","linear model through the origin"))
## You can set the x and y arguments directly to determine where the
## top left corner of the legend is placed on the plot.  Alternately,
## you can set x as "bottomright", "bottom", "bottomleft", "left",
## "topleft", "top", "topright", "right" and "center", and leave y as
## the default (NULL).  More details on how to specify legend
## arguements can be found in the help file.


# use a simple regression model on whale age data
# age and radio isotope ratios from a whale named Moby

TN <- read.table("TeethNitrogen.txt", header = TRUE)
#subset the data to just analyze the data for the single whale "Moby"
head(TN)
summary(TN)
Moby <- subset(TN, Tooth == "Moby")
names(Moby)<- c("N", "Age", "Tooth") 

#plot the data
plot(y = Moby$N, x = Moby$Age,  xlab = "Estimated age", ylab = expression(paste(delta^{15}, "N")) )
#note use of expression() that can be used to include greek symbols

# fit a model to the Moby data set
moby1<- lm(N~Age, data = Moby)
#plot the model
abline(moby1)

#view the model output
moby1

#Evaluate model for: homogeneity, normality, independence
# Let's plot the residuals versus fitted values (homogeneity)
names(moby1)
plot(moby1$fitted, moby1$res)
abline(h=0)

# Histogram and qqnorm plot (normality)
par(mfrow = c(2,1))
hist(moby1$res)
qqnorm(moby1$res)
qqline(moby1$res, distribution = qnorm)

#can also perform Shapiro-Wilk test for normality
shapiro.test(moby1$res)

#residuals versus predictor variable (independence)
par(mfrow = c(1,1))
plot(moby1$res, Moby$Age)


#Diagnostic plots including identification of potential outliers via leverage
plot(moby1)

#evaluate the significance of the regression terms using t-test
summary(moby1)

#evaluate significance of regression terms using F-test
anova(moby1)

plot(y = Moby$N, x = Moby$Age,  xlab = "Estimated age", ylab = expression(paste(delta^{15}, "N")) )
x.seq <- seq(5,50,len=50)

#confidence intervals and predictive intervals
#need to supply newdata to the predict.lm function 
#this must have the same name as the x values in the original model, that is "Age"
(conf.int <- predict.lm(moby1, interval="confidence",
                        newdata=list(Age=x.seq)))
(pred.int <- predict.lm(moby1, interval="prediction",
                        newdata=list(Age=x.seq)))
?predict.lm                             # Note: level=0.95 is the default
## Add to the plot
lines(x.seq, conf.int[,2], lty=2)
lines(x.seq, conf.int[,3], lty=2)
lines(x.seq, pred.int[,2], lty=3)
lines(x.seq, pred.int[,3], lty=3)


#Before moving on, despite seeing little pattern in the residuals vs Age plot,
#what can we say about independence for this analysis? 
# inference to Moby only!

#How do we test for this?
acf(moby1$res)
# temporal autocorrelation which can be modeled with 
#time series models that have an AR(1) or AR(2) structure


#-------   Part II   ----------
## Multiple linear regression

#simulate data
#two covariates and their interaction plus a quadratic of the first covariate 
nobs <- 250
set.seed(123)
x1 <- rnorm(nobs, mean = 0, sd = 1)
x2 <- rnorm(nobs, mean = 0, sd = 1)
#interaction term
x3 <- x1*x2
#quadratic term
x4 <- x1*x1
x0 <- rep(1, times = nobs)

#construct a matrix of covariates
X.mat <- cbind(x0, x1, x2, x3, x4)
head(X.mat)
#vector of regression parameters
beta <- c(5, 2, -1, 0.7, -0.5)

##underlying dynamics
#the quadratic relationship
par(mfrow=c(1,1))
curve(beta[1] + beta[2]*x +beta[5]*x*x, from = -2, to = 2, ylab = "Y", xlab = "X1")

#the interaction effect - a curved surface
x_1 <- seq(-2, 2, length= 30)
x_2<- x_1
f <- function(x_1, x_2) { beta[1] + beta[2]*x_1 + beta[3]*x_2 + beta[4]*x_1*x_2 }
z <- outer(x_1, x_2, f)

op <- par(bg = "white")
persp(x_1, x_2, z, theta = 2, phi = 30, expand = 0.5, col = "lightblue", xlab = "X1", ylab = "X2", zlab = "Mean Obs")

#compare this to no interaction effect - a flat surface
f2 <- function(x_1, x_2) { beta[1] + beta[2]*x_1 + beta[3]*x_2}
z2 <- outer(x_1, x_2, f2)
op <- par(bg = "white")
persp(x_1, x_2, z2, theta = 5, phi = 30, expand = 0.5, col = "pink1", xlab = "X1", ylab = "X2", zlab = "Mean Obs")


#x1 = low
plot(x2, z[,2], type = 'l')

##calculate the mean values using matrix multiplication %*% 
mean.vec <- X.mat %*% beta

##calculate the observed values with measurement error
obs <- mean.vec + rnorm(nobs, mean = 0, sd = 2)

#construct a dataframe to make it easier to model 
sim.df <- data.frame(Obs = obs, X1 = x1, X2 = x2)
pairs(sim.df)

#start modeling with lm
mlm1<-lm(Obs ~ X1 + X2, data = sim.df)

summary(mlm1)

#add terms sequentially to the base model - 
#add an interaction term
mlm2<- lm(Obs ~ X1 + X2 + X1:X2, data = sim.df)
summary(mlm2)

#compare the two models
anova(mlm1, mlm2)

#add a quadratic term for X1, use update()
mlm3 <- update(mlm2, .~. + I(X1^2))
summary(mlm3)

#add a quadratic term for X2
mlm4 <- update(mlm3, .~. + I(X2^2))
summary(mlm4)

#backwards model selection by drop1
drop1(mlm4, test = "F")
#provides an estimate of the SSQ and RSS relative to the base model
#the models are compared the via F test
#AIC computed for each model by default
#suggests that we remove the X2^2 term, which was not one of our simulated covariates

#check of diagnostics
plot(mlm3)

#Evaluation of assumptions
E.sim<- rstandard(mlm3)
#Check for normality
hist(E.sim)
qqnorm(E.sim)
qqline(E.sim)
#Check for independence: residuals versus individual #explanatory variables
plot(y = E.sim, x = sim.df$X1, xlab = "X1", ylab = "Residuals")
abline(0,0)
plot(E.sim, x = sim.df$X2, xlab = "X2", ylab = "Residuals")
abline(0, 0)


#####
## Fit to an actual data set - 
####
#Data set from Appendix A of Zuur et al. 2009
#Forest bird abundances in SE Victoria, Australia
#Covariates are: AREA area of forest patch, DIST distance to patch
#  DISTL distance to larger patch, ALTITUDE altidue of patch, YEAR.ISOL year of clearing, GRAZE index of grazing intensity (levels of 1-5)
Loyn <- read.table("Loyn.txt", header = TRUE)

#create a factor for the GRAZE covariate
Loyn$fGRAZE <- factor(Loyn$GRAZE)
#log transforming several covariates
Loyn$L.AREA<-log10(Loyn$AREA)
Loyn$L.DIST<-log10(Loyn$DIST)
Loyn$L.LDIST<-log10(Loyn$LDIST)

# can use scale() to scale and center covariates as well


Z<-cbind(Loyn$ABUND,Loyn$L.AREA,Loyn$L.DIST,Loyn$L.LDIST,Loyn$YR.ISOL,Loyn$ALT,Loyn$GRAZE)
colnames(Z)<-c("ABUND","L.AREA","L.DIST","L.LDIST","YR.ISOL","ALT","GRAZE")


#define the histogram panel function
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}

#define the panel.cor function
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =  1)
} 

#pairs plot with correlation coefficients 
pairs(Z, upper.panel=panel.cor,diag.panel=panel.hist)

##  Un reto - modelar estos datos con un modelo lineal...
#hint: start with all factors and no interactions

lm1<- lm(ABUND~L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT + fGRAZE, data = Loyn )
summary(lm1)

drop1(lm1, test = "F")

lm2 <- lm(ABUND~L.AREA + fGRAZE, data = Loyn)
summary(lm2)

anova(lm1, lm2)

#to obtain 5 intercepts - 
lm3 <- lm(ABUND~ L.AREA + fGRAZE - 1, data = Lyon)
summary(lm3)





# answers in RetoSolnDay1.R
