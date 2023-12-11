install.packages("zoo")
install.packages("fGarch")
install.packages("stats")
install.packages("base")
install.packages("sm")
install.packages("graphics")
install.packages("xts")
install.packages("ggplot2")
install.packages("sn")
install.packages("quantmod")
install.packages("TTR")
library(quantmod)
library(MASS)
library(fGarch)
library(sn)
library(LambertW)
library(stats)
install.packages("carData")
library(car)
install.packages("knitr")
library(knitr)
install.packages("copula")
library(copula)
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
install.packages("univariateML")
library(univariateML)
install.packages("sm")
library(sm)
install.packages("rugarch")
library(rugarch)
install.packages("parallel")
citation("limma")
library(limma)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages("QRM")
library(QRM)
install.packages("rmarkdown")
library(rmarkdown)
install.packages("moments")
library(moments)
install.packages("MASS")
install.packages("survival")
library(fitdistrplus)
library("evir")
library("forecast")
library("Ecdat")
install.packages("Ecdat")
install.packages("Ecfun")
require("quantmod")

#Download data
#SPXL
getSymbols('SPXL',src='yahoo',from='2012-11-16',to='2022-11-16', periodicity='daily')
summary(SPXL)
#AMZN
getSymbols('AMZN',src='yahoo',from='2012-11-16',to='2022-11-16', periodicity='daily')
summary(AMZN)

#Explore data
#SPXL
#price plot-acf
plot(SPXL$SPXL.Adjusted)
acf(SPXL$SPXL.Adjusted, main="ACF plot of SPXL daily stock prices")#ACF decays to zero slowly. This is a sign of either nonstationarity or possibly of stationarity with long-memory dependence

#Log returns
SPXL.daily.logreturn=na.omit(diff(log(SPXL$SPXL.Adjusted)))
mean(SPXL.daily.logreturn)
skewness(SPXL.daily.logreturn)
kurtosis(SPXL.daily.logreturn)
SPXL.monthly.logreturn=apply.monthly(SPXL.daily.logreturn,Return.cumulative)
acf(SPXL.daily.logreturn, main="ACF plot of SPXL daily logreturn")
# Several of the autocorrelations of the series fall outside the test bounds, which suggests that the series is not white noise.
plot(SPXL.daily.logreturn)
graphics.off()

#Ljung-Box test
Box.test(SPXL.daily.logreturn^2,lag=10, type="Ljung")
# test for white noisehas an extremely small p-value, so the null hypothesis of white noise is strongly rejected. Other choices of K give similar results.
Box.test(SPXL.daily.logreturn,lag=10, type="Ljung")

#AMZN
#price plot-acf
plot(AMZN$AMZN.Adjusted)
acf(AMZN$AMZN.Adjusted, main="ACF plot of AMZN daily stock prices")
#Log returns
AMZN.daily.logreturn=na.omit(diff(log(AMZN$AMZN.Adjusted)))
mean(AMZN.daily.logreturn)
skewness(AMZN.daily.logreturn)
kurtosis(AMZN.daily.logreturn)
plot(AMZN.daily.logreturn)
acf(AMZN.daily.logreturn, main="ACF plot of AMZN daily log returns")
#decays to zero quickly, indicating clearly that the differenced series is stationary.seems no serious correlation but still need to check volativity.good enough?
graphics.off()
#Ljung-Box test
Box.test(AMZN.daily.logreturn^2,lag=5, type="Ljung")#has an extremely small p-value, so the null hypothesis of white noise is strongly rejected. Other choices of K give similar results.
Box.test(AMZN.daily.logreturn^2,lag=510, type="Ljung")

#Time series model
#SPXL
fita=arima(SPXL.daily.logreturn,order=c(1,0,0))
fita
acf(residuals(fita))# not good
acf(abs(residuals(fita)))
SPXL.fit=garchFit(formula=~arma(1,0)+garch(1,1),data=SPXL.daily.logreturn,cond.dist="norm")
summary(SPXL.fit)# alpa1+beta1<1 show it's stationary
coef(SPXL.fit)
SPXL.rsd=residuals(SPXL.fit)#get residual from ar model
SPXL.rsd.std=residuals(SPXL.fit,standardize=TRUE)#get residual from garh model
#Normal distribution
par(mfrow=c(3,2))
acf(SPXL.rsd)
acf(SPXL.rsd^2)
acf(SPXL.rsd.std)
acf(SPXL.rsd.std^2)
qqnorm(SPXL.rsd.std)
qqline(SPXL.rsd.std,col="red")

#t distribution
td.SPXL=fitdistr(SPXL.rsd.std,"t")
qqt(SPXL.rsd.std,df=td.SPXL$estimate[3],ylim=range(SPXL.rsd.std),main = "t Q-Q Plot, SPXL", 
   xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", plot.it = TRUE)
qqline(SPXL.rsd.std,col="red")
graphics.off()

#AMZN
fita=arima(AMZN.daily.logreturn,order=c(1,0,0))
fita
acf(residuals(fita))# seem good
acf(abs(residuals(fita)))
AMZN.fit=garchFit(formula=~arma(1,0)+garch(1,1),data=AMZN.daily.logreturn,cond.dist="norm")
coef(AMZN.fit)
AMZN.rsd=residuals(AMZN.fit)
AMZN.rsd.std=residuals(AMZN.fit,standardize=TRUE)
#Normal distribution
par(mfrow=c(3,2))
acf(AMZN.rsd)
acf(AMZN.rsd^2)
acf(AMZN.rsd.std)
acf(AMZN.rsd.std^2)
qqnorm(AMZN.rsd.std)
qqline(AMZN.rsd.std,col="red")

#t distribution
td.AMZN=fitdistr(AMZN.rsd.std,"t")
qqt(AMZN.rsd.std,df=td.AMZN$estimate[3],ylim=range(AMZN.rsd.std),main = "t Q-Q Plot, AMZN", 
    xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", plot.it = TRUE)
qqline(AMZN.rsd.std,col="red")
graphics.off()

#Copulas
#SPXL
u.fit.SPXL=as.numeric(pt(SPXL.rsd.std, df=td.SPXL$estimate[3]))
#AMZN
u.fit.AMZN=as.numeric(pt(AMZN.rsd.std, df=td.AMZN$estimate[3]))
#Kendall's tau
tau=cor(u.fit.SPXL, u.fit.AMZN, method = "kendall")
tau
#Spearman's Rho
omega=sin(tau*pi/2)
omega

U.hat=cbind(u.fit.SPXL,u.fit.AMZN)

#t-copula
Ct=fitCopula(copula=tCopula(dim=2),data=U.hat,method="ml",start=c(omega,6));#fit t copula
Ct@estimate;
loglikCopula(param=Ct@estimate,u=U.hat,copula=tCopula(dim=2));#compute loglikelihood function
-2*.Last.value+2*length(Ct@estimate);#compute AIC

#Gaussian copula
Cgauss=fitCopula(copula=normalCopula(dim=2),data=U.hat,method="ml",start=c(omega));#fit Gaussian copula
Cgauss@estimate;
loglikCopula(param=Cgauss@estimate,u=U.hat,copula=normalCopula(dim=2));
-2*.Last.value+2*length(Cgauss@estimate);#compute AIC

#Frank copula
Cfr=fitCopula(copula=frankCopula(1,dim=2),data=U.hat,method="ml");#fit frank copula
Cfr@estimate;
loglikCopula(param=Cfr@estimate,u=U.hat,copula=frankCopula(dim=2));
-2*.Last.value+2*length(Cfr@estimate);#compute AIC

#Clayton copula
Ccl = fitCopula(copula=claytonCopula(1, dim=2), data=U.hat, method="ml")
Ccl@estimate
loglikCopula(param=Ccl@estimate, u=U.hat, copula=claytonCopula(dim = 2))
-2*.Last.value + 2*length(Ccl@estimate)

#Gumbel copula
Cgu=fitCopula(copula=gumbelCopula(dim=2),data=U.hat,method="ml")
Cgu@estimate
loglikCopula(param=Cgu@estimate,u=U.hat,copula=gumbelCopula(dim=2))
#AIC
-2*.Last.value+2*length(Cgu@estimate)

#Sample
n=1000
Simu_U=rCopula(n,tCopula(dim=2,Ct@estimate[1],df=Ct@estimate[2]))
Simu_U
plot(Simu_U[,1],Simu_U[,2])
#transform marginals into the estimated t
Simu_X1=qt(Simu_U[,1],df=td.SPXL[["estimate"]][["df"]])
Simu_X2=qt(Simu_U[,2],df=td.AMZN[["estimate"]][["df"]])
plot(Simu_X1,Simu_X2)
Simu_X=cbind(Simu_X1,Simu_X2) #simulated sample from the fitted marginals and copula

#Risk calculation
mu1=SPXL.fit@fit[["par"]][["mu"]]
mu2=AMZN.fit@fit[["par"]][["mu"]]
sigma1=sqrt(SPXL.fit@fit[["par"]][["omega"]]+SPXL.fit@fit[["par"]][["alpha1"]]*SPXL.rsd[length(SPXL.rsd)]^2+SPXL.fit@fit[["par"]][["beta1"]]*SPXL.fit@sigma.t[length(SPXL.fit@sigma.t)]^2)
sigma2=sqrt(AMZN.fit@fit[["par"]][["omega"]]+AMZN.fit@fit[["par"]][["alpha1"]]*AMZN.rsd[length(AMZN.rsd)]^2+AMZN.fit@fit[["par"]][["beta1"]]*AMZN.fit@sigma.t[length(AMZN.fit@sigma.t)]^2)
xi1=SPXL.fit@fit[["par"]][["ar1"]]
xi2=AMZN.fit@fit[["par"]][["ar1"]]
for(i in c(1:9)/10) {#consider rho=0.1,0.2...,0.9
  VaR_=(i*(mu1+xi1*as.numeric(SPXL.daily.logreturn[length(SPXL.daily.logreturn)])+sigma1*Simu_X1))+((1-i)*(mu2+xi2*as.numeric(AMZN.daily.logreturn[length(AMZN.daily.logreturn)])+sigma2*Simu_X2))
  VaR99=quantile(VaR_,0.99)
  cat("rho=",i,",VaR=",VaR99,"\n")
  }
  
