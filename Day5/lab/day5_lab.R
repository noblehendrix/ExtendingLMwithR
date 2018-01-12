### ------------------------------------------------------------
## QERM 514 Lab 9 | Spring 2014 | University of Washington
## Cole Monnahan | monnahc@uw.edu | 5/30/2014

## Generalized Additive Models (GAMs)

##  Objectives:
##	- Learn how to build an additive model in R
##	- Compare GAM models
##	- Use various functions associated with GAM models
##      - Use a smoother to inform a possible parametric form
### ------------------------------------------------------------


## Places to read about additive models:
### Faraway ELM Ch. 12
### Generalized Additive Models: An Introduction With R
### Ch. 34 (Crawley) linked on course website

## Read in the data
library(faraway)
data(ozone)
head(ozone)

## O3 = atmospheric ozone concentration
## temp = temperature measured at El Monte (suburb, east LA)
## ibh = inversion base height at LAX
#	(mixing height, estimated height where pollutants from surface mix with air)
## ibt = inversion top temperature at LAX
#	http://en.wikipedia.org/wiki/Inversion_%28meteorology%29

## Have a look at the data so we know what we're working with. Notice
## that you can pass different functions to different parts of the
## panels. Here I'm specifying to add a loess curve (panel.smooth) in
## the lower half, and the correlations in the upper half. Note that
## histograms on the diagonal can be used with diag.panel=panel.hist,
## but this precludes specifying your own color (unless you edit the
## function!). Note the use of transparency using the
## rgb(red,blue,green, alpha) function.

## We're only going to look at a few of these:
pairs(O3~ temp+ ibh + ibt, data=ozone, col=rgb(0,0,0,.5), pch=16)

## Start with a linear model (a base to compare the GAM models to)
mod.lm <- lm(O3 ~ temp + ibh + ibt, data=ozone)
summary(mod.lm)

## Now let's look at the GAM versions of the model

## We will look at two packages that allow fitting via the gam()
## function (1) gam and (2) mgcv.

## other options include: gss, assist, gamlss, vgam (for more advanced modeling)
## install.packages("gam")
detach("package:mgcv", unload=TRUE)
library(gam)

## GAM model with loess smoother
mod.gam1 <- gam(O3 ~ lo(temp) + lo(ibh) + lo(ibt), data=ozone)

## Faraway: "The p-values are only approximate at best and should be
## viewed with some skepticism. It is generally better to fit the
## model without the predictor of interest and then construct the
## F-test" (see below for this)
summary(mod.gam1)

## Test the significance of ibt for the model:
mod.gam2 <- gam(O3 ~ lo(temp) + lo(ibh), data=ozone)
anova(mod.gam2, mod.gam1, test="F") #p-value still approximate
## F-test says do not reject the null that lo(ibt) = 0, i.e. ibt is
## not significant

mod.gam3 <- gam(O3 ~ lo(ibh) + lo(ibt), data=ozone)
anova(mod.gam1, mod.gam3, test="F")
## Temp is highly significant! This may be more apparent when looking
## at the model fit

par(mfrow=c(1,3))
plot.gam(x=mod.gam1, residuals=TRUE, se=TRUE, pch=16, col=rgb(0,0,0,.5),
         cex=.5, scale=40)
## Residuals = TRUE adds the points, se=TRUE adds the dotted lines
## (set rug=F if you don't want the tick marks on the bottom)
## from plot:  temp changes slope at 60, ibh has max at 1000, ibt may be constant

## From looking at temp it might be that a simple linear predictor
## would be just as good.
mod.gam4 <- gam(O3 ~ temp + lo(ibh), data=ozone)
AIC(mod.gam4, mod.gam2)
## Apparently not!

## We can compare AIC values between GAM and simple linear models:
AIC(mod.lm, mod.gam1, mod.gam2, mod.gam3, mod.gam4)

## Some of the things you can extract from a gam object:
mod.gam1$coeff #or coeff(mod.gam1)
## "the coefficients of the parametric part of the additive predictors,
## which multiply the columns of the model matrix."
fit.add <- mod.gam1$additive
##"the additive fit, given by the product of the model matrix and the
## coefficients, plus the columns of the smooth component."
fitted <- mod.gam1$fitted #or predict(mod.gam1)
##"the fitted mean values, obtained by transforming the component
## additive.predictors using the inverse link function"

## Let's plot these vs the observed values to get a feel for what is
## going on:
par(mfrow=c(1,2))
plot(ozone$O3,fit.add)
plot(ozone$O3, fitted)
## They are the same because we are using the default Gaussian family
## (i.e. this is an additive model not a GAM!).

gam.resid <- mod.gam1$residuals #"typically not interpretable without
                                        #rescaling the weightsl"
lm.resid <- residuals(mod.lm)
plot(ozone$O3, gam.resid)
plot(ozone$O3, lm.resid, col="red")

## Look up gam package for more!
?gam

## Calculate R^2 for this model:
summary(mod.gam1)
1 - 5935.096 / 21115.41


## Sidenote: A toy example built into gam.
data(gam.data)
pairs(y~ x +z, gam.data)
## Choose the "target equivalent df" of the smoother:
gam.object <- gam(y ~ s(x,6) + z,data=gam.data)
summary(gam.object)
par(mfrow=c(1,2))
plot(gam.object,se=TRUE, residuals=TRUE, scale=5)
data(gam.newdata); gam.newdata
predict(gam.object, newdata=gam.newdata)
gam.object2 <- gam(y ~ s(x,6),data=gam.data)
anova(gam.object2, gam.object, test="F")
AIC(gam.object2, gam.object)
## ------------------------------------------------------------

## Refit the model using another package
## install.packages("mgcv")
?gam
detach("package:gam")                   # detach gam to be safe
?gam
library(mgcv)
?gam
## mgcv only allows splines for the smoother and the appropriate level
## of smoothing is chosen by default

## more restrictive = avoids subjectivity, negative= auto select can fail

## This now gives an error since mgcv won't let us pick a smoother df.
gam.object <- gam(y ~ s(x,6) + z,data=gam.data)
## So refit and examine summary
gam.object <- gam(y ~ s(x) + z,data=gam.data)
summary(gam.object)                     # Notice the difference!

## Fit our GAM:
mod.mgcv1 <- gam(O3 ~ s(temp) + s(ibh) + s(ibt), data=ozone)
summary(mod.mgcv1) # ibt doesn't appear to be significant
mod.mgcv1$aic

## Examine the fit:
par(mfrow=c(1,3))
plot(mod.mgcv1, residuals=TRUE, pch=".",  lwd=1.5)
## shows the same pattern as above (Note: "partial residuals" =
## Pearson residuals added to the smooth terms evaluated at the
## appropriate covariate values. In a well-fitting model, the partial
## residuals should be evenly scattered around the curve.)

## Is there really is a change in the trend for temp?
am1 <- gam(O3 ~ s(temp) + s(ibh), data=ozone)
am2 <- gam(O3 ~ temp + s(ibh), data=ozone)
anova(am2, am1, test="F") #p-value still approx. but looks like Yes, smoother is sig.
## (Note: also significant when we include s(ibt) in the models)
## (Note: You can also do anova(am1, mod.mgcv1, test="F"), borderline NS)

## Predict ozone levels from the model, with Standard Errors
predict(mod.mgcv1, newdata=data.frame(temp=60, ibh=2000, ibt=100), se=T)

## Diagnostic plots
gam.check(mod.mgcv1)
par(mfrow=c(1,1))
plot(predict(mod.mgcv1), residuals(mod.mgcv1)); abline(h=0)
## Nonconstant variance...

## Just for kicks, let's re-run it as poisson
pois.mgcv1 <- gam(O3 ~ s(temp) + s(ibh) + s(ibt), family=poisson,
                  scale=-1, data=ozone)
summary(pois.mgcv1) #negative scale indicates scale parameter is unknown

## Check what the smoothing functions look like under the poisson model
par(mfrow=c(1,3))
plot(pois.mgcv1, residuals=TRUE)

pois.mgcv2 <- gam(O3 ~ s(temp) + s(ibh), family=poisson, scale=-1,
                  data=ozone)
anova(pois.mgcv2, pois.mgcv1, test="F") #NS so we could take out ibt
AIC(pois.mgcv1, pois.mgcv2) # not very different
gam.check(pois.mgcv2)

## How does the AM vs GAM look?
par(mfrow=c(1,2))
plot(predict(mod.mgcv1), residuals(mod.mgcv1)); abline(h=0)
plot(predict(pois.mgcv2), residuals(pois.mgcv2)); abline(h=0) #looks better


## Also note that it is possible to fit *interaction* effects with a
## gam as follows. See Faraway ELM p.236 for a discussion on why this
## model is not appropriate. It is hard to tell from the summary and
## F-tests, so actually do it!
mod.mgcvint <- gam(O3~s(temp,ibh)+s(ibt), data=ozone)
plot(mod.mgcvint)

## Another thing you can play around with (see help files):
vis.gam(pois.mgcv2)

### ------------------------------------------------------------
## Using a smoother to guide our model selection process for parametric
## terms

## One use of the GAM is as an exploratory tool.  Use the GAM plots to add
## "broken stick" regression to our original linear model
par(mfrow=c(1,3))
plot(mod.mgcv1, residuals=TRUE)
## Our GAM plots suggest a change in slope for temp at 60 degrees
## And a change in slope for ibh at 1000

## One way to do this:
## Create new variables that begin their count at the "break point"
temp2 <- (ozone$temp>60)*(ozone$temp-60)
## Create a new vector that only counts temps above 60, and subtract
## 60 to begin count at 1, 2, 3...

temp2			# Look at what you've created
## and compare to the temp variable
par(mfrow=c(1,1))
plot(ozone$temp, temp2)

##  Do the same thing for ibh
ibh2 <- (ozone$ibh-1000) * (ozone$ibh>1000)

## Now run the "broken stick" regression
mod.lm2 <- lm(O3 ~ temp + temp2 + ibh + ibh2, data=ozone)
summary(mod.lm2)

## So what have we learned from the GAM? does the broken stick work
## better?
mod.lm3 <- lm(O3 ~ temp + ibh, data=ozone)
AIC(mod.lm2, mod.lm3, mod.mgcv1)

## Yes! Would you have tried this without fitting a GAM? Probably
## not. In this way the GAM has taught us something about the
## underlying structure of the model.

## End of lab
### ------------------------------------------------------------

