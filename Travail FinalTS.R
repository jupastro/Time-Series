
#Load the required packages:
require(forecast)
library(tibbletime)
library(dplyr)
require(xts)
require(TSstudio)
require(fGarch)
library(rugarch)
library(tseries)
library(fBasics)
require(aTSA)
source(cd) ## needs to be in the working directory
#We create the time series 
Stock <- ts(read.table('dell.txt', header = F)$V1); ## or "import Dataset" in Rstudio
basicStats(Stock)# We obtain an aproximation of the mean and the variance

# Stock
# nobs        1261.000000
# NAs            0.000000
# Minimum      -15.741931
# Maximum       18.849372
# 1. Quartile   -1.823277
# 3. Quartile    2.482740
# Mean           0.365698
# Median         0.365047
# Sum          461.144737
# SE Mean        0.094782
# LCL Mean       0.179749
# UCL Mean       0.551646
# Variance      11.328433
# Stdev          3.365774
# Skewness      -0.052559
# Kurtosis       1.715090  we can see that the distribution would be more centered than an standard normal one


shapiro.test(Stock)
# Shapiro-Wilk normality test
# 
# data:  Stock
# W = 0.98584, p-value = 9.801e-10 as the p value of the test is <<< than p=0.05 thus we night reject  the null hypothesis that the series is based on a normal distribution 

daily <- (ts(Stock, start = c(1993,304*253/356), frequency = 253))
# we create a time series starting from 1993 and with trading year frequency
# We use 253 as they are the trading days per year
daily.squared=daily^2;
# As we are dealing with highly volatile finantial data it is important to take a look to the squared ts
#Plot both:
par(mfrow=c(1,1))
plot(daily,main="Daily stock returns ") ## plot time series which is actually a previous log
plot(daily.squared,main='Square of the daily stock returns')
#Now we must check if any kind of ARMA model fits for the daily stock returns

#Check if any of them has a mean tendence:
par(mfrow=c(1,1))
plot(daily,main="Daily stock returns ") ## plot time series
t.daily <- seq(1993 + 304/365, by = 1/253, length = 1261)
abline(reg=lm(daily~t.daily),col=2)

#Check for stationary conditions
par(mfrow=c(2,1))
plot(daily[1:253], type="l", ylim=c(min(daily), max(daily)),main="Comparison of the daily stock returns")
# we plot every year one compared to each other so that we can establish their trends
for(i in 2:5){
  lines(daily[253*(i-1)+(1:253)],col=i)
}
legend("bottomright", legend = 1994:1999, lty = 1, col = 1:12, cex = 0.5)

plot(daily.squared[1:253], type="l", ylim=c(min(daily.squared), max(daily.squared)),main="Comparison of the daily squared stock returns")
# we plot every year one compared to each other so that we can establish their trends
for(i in 2:5){
  lines(daily.squared[253*(i-1)+(1:253)],col=i)
}
legend("bottomright", legend = 1994:1999, lty = 1, col = 1:12, cex = 0.5)

#Check correlations for a white noise:
?acf
par(mfrow=c(2,1))
acf=acf(daily, lag = 50, main = 'ACF of daily stock returns')
par(new=TRUE)
plot(acf,ci=0.99)

pacf=pacf(daily, lag = 50, main = 'PACF of daily stock returns')
par(new=TRUE)
plot(pacf,ci=0.99)
# Looks like there's some partial autocorrelation for the lag 2
# Fit an ARIMA model:
Model_Select=expand.grid(order_p=0:4,order_q =0:4,Aic=0,Bic=0)
n=1
while(n<length(Model_Select$order_p)+1){
  suppressWarnings( model <- arima(daily,order=c(Model_Select$order_p[n],0,Model_Select$order_q[n])))
  Model_Select$Aic[n] <- Aic.arima(daily, model ) ## AIC order i
  Model_Select$Bic[n] <- Bic.arima(daily, model) ## BIC order i
  n=n+1
}
suppressWarnings(Comp.Sarima(daily, d=0, saison=0, D=0, p.max=4, q.max=4, P.max=0, Q.max=0))
Best_Aic=which.min(Model_Select$Aic)
Best_Bic=which.min(Model_Select$Bic)
modelAIC <- arima(daily,order=c(Model_Select$order_p[Best_Aic] ,0,Model_Select$order_q[Best_Aic] ))
modelBIC <- arima(daily,order=c(Model_Select$order_p[Best_Bic] ,0,Model_Select$order_q[Best_Bic] )) 
modelAr1<- arima(daily,order=c(1,0,0))
#Check the importance of the coefficients 
coef.BIC<-coef.p(modelBIC$coef, diag(modelBIC$var.coef))
coef.AIC<-coef.p(modelAIC$coef, diag(modelAIC$var.coef))
coef.ar1<-coef.p(modelAr1$coef,diag(modelAr1$var.coef))
coef.BIC
coef.AIC
coef.ar1
# Check the residuals as white noise
tsdiag(modelBIC,main='BIC')
tsdiag(modelAIC,main='AIC')
par(mfrow=c(2,1))
acf(modelAIC$res, lag = 50, main = 'ACF of residuals model AIC')
pacf(modelAIC$res, lag = 50, main = 'PACF of residuals model AIC')
par(mfrow=c(2,1))
acf(modelBIC$res, lag = 50, main = 'ACF of residuals model BIC')
pacf(modelBIC$res, lag = 50, main = 'PACF of residuals model BIC')

qqnorm(modelAIC$res) ## check for normality of residuals with a QQplot
abline(0, 23/7, col = 2)
qqnorm(modelBIC$res) ## check for normality of residuals with a QQplot
abline(0, 23/7, col = 2)

?tsdiag()
#Check the squared residuals as white noise
par(mfrow=c(2,1))
acf((modelAIC$res)^2, lag = 50, main = 'ACF of squared residuals model AIC')
pacf((modelAIC$res)^2, lag = 50, main = 'PACF of  squared residuals model AIC')
par(mfrow=c(2,1))
acf((modelBIC$res)^2, lag = 50, main = 'ACF of squared residuals model BIC')
pacf((modelBIC$res)^2, lag = 50, main = 'PACF of  squared residuals model BIC')
par(mfrow=c(1,2))


#Fit a (G)ARCH model to our two previous models:
suppressWarnings(Comp.Sarima(modelAIC$residuals^2, d=0, saison=0, D=0, p.max=4, q.max=4, P.max=0, Q.max=0))
# mod?le (p,d,q)x(P,D,Q)_saison :  3 0 3 x 0 0 0 _ 0 :  nb param:  6     AIC: 0 
# mod?le (p,d,q)x(P,D,Q)_saison :  3 0 4 x 0 0 0 _ 0 :  nb param:  7     AIC: 10.59227 
# mod?le (p,d,q)x(P,D,Q)_saison :  4 0 4 x 0 0 0 _ 0 :  nb param:  8     AIC: 3.783044 
suppressWarnings(Comp.Sarima(modelBIC$residuals^2, d=0, saison=0, D=0, p.max=4, q.max=4, P.max=0, Q.max=0))
# mod?le (p,d,q)x(P,D,Q)_saison :  2 0 4 x 0 0 0 _ 0 :  nb param:  6     AIC: 12.67532 
# mod?le (p,d,q)x(P,D,Q)_saison :  3 0 3 x 0 0 0 _ 0 :  nb param:  6     AIC: 0 
# mod?le (p,d,q)x(P,D,Q)_saison :  4 0 4 x 0 0 0 _ 0 :  nb param:  8     AIC: 5.767767 
#Thus we will use an ARMA(3,3) to model the residuals which implies and GARCH(1-3,3)
n=1
#With a simple normal distribution 
ARMA.GARCH.1 <- garchFit(data ~ arma(2,0) + garch(1,1), data = daily, trace = F)# We have reduced the order due to the analysis of the coefficients 
capture.output(summary(ARMA.GARCH.1)) ## print only the interesting part// best one according to AIC and BIC 
ARMA.GARCH.2 <- garchFit(data ~ arma(2,0) + garch(1,2), data = daily, trace = F)
 capture.output(summary(ARMA.GARCH.2)) ## print only the interesting part
ARMA.GARCH.3 <- garchFit(data ~ arma(2,0) + garch(1,3), data = daily, trace = F)
 capture.output(summary(ARMA.GARCH.3)) ## print only the interesting part

 ARMA.GARCH.4 <- garchFit(data ~ arma(0,0) + garch(1,1), data = daily, trace = F)
 capture.output(summary(ARMA.GARCH.4)) ## print only the interesting part
 ARMA.GARCH.5 <- garchFit(data ~ arma(0,0) + garch(1,2), data = daily, trace = F)
 capture.output(summary(ARMA.GARCH.5)) ## print only the interesting part
 ARMA.GARCH.6 <- garchFit(data ~ arma(0,0) + garch(1,3), data = daily, trace = F)
 capture.output(summary(ARMA.GARCH.6)) ## print only the interesting part
 
 #We apply a conditional t-student distribution
  ?garchFit
 ARMA.GARCH.1 <- garchFit(data ~ arma(2,0) + garch(1,1), data = daily, trace = F,include.skew = FALSE,include.shape =FALSE,cond.dist='sstd')# We have reduced the order due to the analysis of the coefficients 
 capture.output(summary(ARMA.GARCH.1)) ## print only the interesting part// best one according to AIC and BIC 
 ARMA.GARCH.2 <- garchFit(data ~ arma(2,0) + garch(1,2), data = daily, trace = F,include.skew = FALSE,include.shape =FALSE,cond.dist='sstd')
 capture.output(summary(ARMA.GARCH.2)) ## print only the interesting part
 ARMA.GARCH.3 <- garchFit(data ~ arma(2,0) + garch(1,3), data = daily, trace = F,include.skew = FALSE,include.shape =FALSE,cond.dist='sstd')
 capture.output(summary(ARMA.GARCH.3)) ## print only the interesting part
 
 ARMA.GARCH.4 <- garchFit(data ~ arma(0,0) + garch(1,1), data = daily, trace = F,include.skew = FALSE,include.shape =FALSE,cond.dist='sstd')
 capture.output(summary(ARMA.GARCH.4)) ## print only the interesting part
 ARMA.GARCH.5 <- garchFit(data ~ arma(0,0) + garch(1,2), data = daily, trace = F,include.skew = FALSE,include.shape =FALSE,cond.dist='sstd')
 capture.output(summary(ARMA.GARCH.5)) ## print only the interesting part
 ARMA.GARCH.6 <- garchFit(data ~ arma(0,0) + garch(1,3), data = daily, trace = F,include.skew = FALSE,include.shape =FALSE,cond.dist='sstd')
 capture.output(summary(ARMA.GARCH.6)) ## print only the interesting part
 

 plot((abs(daily)))
 1#Interesting till here # Best model AIC-log-likelihood: ARMA.GARCH1; if we want the most parsimonious one GARCH(1,1)
 # Now we must check the residuals 
par(mfrow=c(2,1))
acf(ARMA.GARCH.1@residuals,main='AR(2)xGARCH(1,1) residuals',ci=0.99)
pacf(ARMA.GARCH.1@residuals,main='AR(2)xGARCH(1,1) residuals',ci=0.99)

par(mfrow=c(2,1))
acf(ARMA.GARCH.4@residuals, main='GARCH(1,1) residuals',ci=0.99)
pacf(ARMA.GARCH.4@residuals,main='GARCH(1,1) residuals',ci=0.99)

# now we check for the standardized residuals

par(mfrow=c(3,5))
acf(ARMA.GARCH.1@residuals^2/ARMA.GARCH.1@h.t,main='AR(2)xGARCH(1,1) standardized squared residuals',ci=0.99)
pacf(ARMA.GARCH.1@residuals^2/ARMA.GARCH.1@h.t,main='AR(2)xGARCH(1,1) standardized squared residuals',ci=0.99)

plot(ARMA.GARCH.1,which='all')


acf(ARMA.GARCH.4@residuals^2/ARMA.GARCH.1@h.t, main='GARCH(1,1) standardized squared residuals',ci=0.99)
pacf(ARMA.GARCH.4@residuals^2/ARMA.GARCH.1@h.t,main='GARCH(1,1) standardized squared residuals',ci=0.99)

 #Comparison of their prediction behaviour 

par(mfrow=c(2,1))
model<-ugarchspec(variance.model = list(model = "fGARCH",submodel='GARCH', garchOrder = c(1, 1)), 
                  mean.model = list(armaOrder = c(2, 0), include.mean = TRUE), distribution.model = "sstd")
modelfit<-ugarchfit(spec=model,data=daily)

spec = getspec(modelfit);
setfixed(spec) <- as.list(coef(modelfit));
length(daily)
forecast = ugarchforecast(spec, n.ahead = 1, n.roll = 1260, data = daily[1:1261], out.sample = 1260);

plot(forecast,which=4)
model<-ugarchspec(variance.model = list(model = "fGARCH",submodel='GARCH', garchOrder = c(1, 1)), 
                  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "sstd")
modelfit<-ugarchfit(spec=model,data=daily)

spec = getspec(modelfit);
setfixed(spec) <- as.list(coef(modelfit));
length(daily)
forecast = ugarchforecast(spec, n.ahead = 1, n.roll = 1260, data = daily[1:1261,drop=FALSE], out.sample = 1260);

plot(forecast,which=4)

#Future prediction:
ARMA.GARCH.1 <- garchFit(data ~ arma(2,0) + garch(1,1), data = daily[1:round(1*length(daily))], trace = F)# We have reduced the order due to the analysis of the coefficients 
par(mfrow=c(2,2))
pred1=fGarch::predict(ARMA.GARCH.1,n.ahead=round(0.1*length(daily)),plot=T,nx=70)
ARMA.GARCH.1 <- garchFit(data ~ arma(0,0) + garch(1,1), data = daily[1:round(1*length(daily))], trace = F)# We have reduced the order due to the analysis of the coefficients 
pred2=fGarch::predict(ARMA.GARCH.1,n.ahead=round(0.1*length(daily)),nx=70,plot=T)

ARMA.GARCH.1 <- garchFit(data ~ arma(2,0) + garch(1,1), data = daily[1:round(0.98*length(daily))], trace = F)# We have reduced the order due to the analysis of the coefficients 

pred1=fGarch::predict(ARMA.GARCH.1,n.ahead=round(0.1*length(daily)),plot=T,nx=70)
ARMA.GARCH.1 <- garchFit(data ~ arma(0,0) + garch(1,1), data = daily[1:round(0.98*length(daily))], trace = F)# We have reduced the order due to the analysis of the coefficients 
pred2=fGarch::predict(ARMA.GARCH.1,n.ahead=round(0.1*length(daily)),nx=70,plot=T)



#Volatility prediction

par(mfrow=c(2,1))
model<-ugarchspec(variance.model = list(model = "fGARCH",submodel='GARCH', garchOrder = c(1, 1)), 
                  mean.model = list(armaOrder = c(2, 0), include.mean = TRUE), distribution.model = "norm")
modelfit<-ugarchfit(spec=model,data=daily)

spec = getspec(modelfit);
setfixed(spec) <- as.list(coef(modelfit));
length(daily)
forecast = ugarchforecast(spec, n.ahead = round(0.2*length(daily)),data = daily[1:1261]);

plot(forecast,which=3)
model<-ugarchspec(variance.model = list(model = "fGARCH",submodel='GARCH', garchOrder = c(1, 1)), 
                  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "norm")
modelfit<-ugarchfit(spec=model,data=daily)

spec = getspec(modelfit);
setfixed(spec) <- as.list(coef(modelfit));
length(daily)
forecast = ugarchforecast(spec, n.ahead = round(0.2*length(daily)), data = daily[1:1261,drop=FALSE]);

plot(forecast,which=3)

 
