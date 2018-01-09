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

#start with all factors and no interactions
M1 <- lm(ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT +
           fGRAZE, data = Loyn)
summary(M1)

anova(M1)

#drop the fGRAZE factor
M2 <- lm(ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT, data = Loyn)
anova(M1, M2)

#use a backwards selection by dropping 1 factor at a time
drop1(M1,test="F")

#a selected model given the backwards selection approach
M3 <- lm(ABUND ~ L.AREA + fGRAZE, data = Loyn)

#run an automatic stepwise backwards selection using AIC
step(M1)

#comparison of the reduced model with the full model
anova(M1, M3)  

plot(M3) #standard graphical output
op <- par(mfrow = c(2, 2))
#Check for normality
E <- rstandard(M3)
hist(E)
#qqnorm(E)
#Check for independence: residuals versus individual #explanatory variables
plot(y = E, x = Loyn$L.AREA, xlab = "AREA", ylab = "Residuals")
abline(0,0)
plot(E ~ Loyn$fGRAZE, xlab = "GRAZE", ylab = "Residuals")
abline(0, 0)
par(op)

plot(Loyn$L.AREA,Loyn$ABUND)
#Predictions - 
#make some predictions from the model: 
#Loyn$SL.AREA<-sort(Loyn$L.AREA)
D1<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==1],fGRAZE="1")
D2<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==2],fGRAZE="2")
D3<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==3],fGRAZE="3")
D4<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==4],fGRAZE="4")
D5<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==5],fGRAZE="5")

P1<-predict(M3,newdata=D1)
P2<-predict(M3,newdata=D2)
P3<-predict(M3,newdata=D3)
P4<-predict(M3,newdata=D4)
P5<-predict(M3,newdata=D5)


D1<-data.frame(L.AREA = Loyn$L.AREA, fGRAZE = "1")
P1<-predict(M3, newdata = D1)
lines(D1$L.AREA,P1,lty=1)
lines(D2$L.AREA,P2,lty=2)
lines(D3$L.AREA,P3,lty=3)
lines(D4$L.AREA,P4,lty=4)
lines(D5$L.AREA,P5,lty=5)
