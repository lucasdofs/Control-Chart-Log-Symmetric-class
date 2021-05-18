############### Author: LUCAS SALES ############
#Institution: Universidade Federal do Rio Grande do Norte
#e-mail: lucasofsales@gmail.com  or ldo-sales@hotmail.com
#Log-symmetrical Class simulation #######


#This code aims to present a monitoring of the PVI of engine oils(see the attached csv set)
#using the method proposed by Lucas Sales, André Pinho,Francisco Medeiros and 
#Marcelo Bourguignon (paper submitted to the Statistical Methods & Applications)
library("gamlss")
library(psych)
library(stats)
library(ICSNP)
library(xtable)
library(LaplacesDemon) 
library(e1071)
library(astsa)

oil=read.table("oil_quality.csv",sep=",")
teste=c(oil$V1,oil$V2,oil$V3)
oil=as.matrix(oil)
mean(teste)
sd(teste)
hist(teste,main="",xlab= "PVIS",xlim=c(60,160))
boxplot(teste,main="", ylab= "PVIS")

#Verifying if the Lognormal really fitted best for the data than the usual distribution
#for non-negativo assymetric data
fit1=gamlss(teste~1, family=LOGNO)
fit2=gamlss(teste~1, family=WEI())
fit3=gamlss(teste~1, family=EXP())
fit4=gamlss(teste~1, family=GA())

fit1$aic;fit2$aic;fit3$aic;fit4$aic

#Checking the residuals
qqPlot(fit1$residuals,envelope=.99) #Normality ok 
acf2(fit1$residuals) #No autocorrelation

med.hat=4.589
sig.hat=0.167

media.hat=mean(teste)
sd.hat=sd(teste)/sqrt(3)

PEI=function(x){
  sqrt(mean(x)*harmonic.mean(x))
}

med.oil1=NULL
for(i in 1:25){
  med.oil1[i]=PEI(oil[i,1:3])
}
med.oil1

#Obtainig control limits

amostraPEI=function(n,med,phi,delta){
  x=NULL
  for(i in 1:10000){
    # print(i)
    x[i]=PEI(rlnorm(n,meanlog=log(med*(1+delta)),sdlog=phi))
  }
  return(x)
}

amostraHL=function(n,med,phi,d,delta){
  x=NULL
  for(i in 1:10000){
    # print(i)
    x[i]=hl.loc(rlnorm(n,meanlog=log(med*(1+delta)),sdlog=phi))/(1+d)
  }
  return(x)
}

alpha=0.00198
x0=amostraPEI(3,exp(med.hat),sig.hat,0)
z0=(x0-exp(med.hat))/sd(x0)
(q1=quantile(z0,probs=c(alpha/2)))
(q2=quantile(z0,probs=c(1-alpha/2)))
fda0=ecdf(z0)
#ARL0 
(arl0=1/(fda0(q1)+1-fda0(q2))) #500
(arl0f=1/(fda0(-3.09)+1-fda0(3.09))) #370.35


LIC=exp(med.hat)+q1*sd(med.oil1) 
LSC=exp(med.hat)+q2*sd(med.oil1) 

ULIC=98.47-3.0*sd(med.oil1) 
ULSC=98.47+3.0*sd(med.oil1)

delta=0.00
#x1=amostraPEI(3,exp(med.hat),sig.hat,delta)
#rlnorm(30,meanlog=(med.hat*(1+delta)),sdlog=sig.hat)

#x1 is the random samples of phase 2 using in the original paper
x1=c(99.06930, 113.93315,  93.44010, 109.76233, 109.69995,  84.88708,
     100.94031, 111.77178,  90.35337,  86.78346,  89.70118, 101.73825,
     129.61654, 92.41196, 100.68812, 102.76927, 106.72422, 105.36035,
     92.16729, 105.04633)
plot(c(med.oil1,x1[-13]),pch=16,col="blue",ylim=c(65,170),xlab="Sample",
     ylab="Median of the PVIS" )
#points(x=38,y=x1[13],pch=16,col="red")
points(x=38,y=x1[13],pch=8,col="red")

abline(h=LIC,col="red")
abline(h=LSC+0.1,col="red")

abline(h=ULIC,col="blue",lty=2)
abline(h=ULSC,col="blue",lty=2)

lines(x=c(25,25),y=c(ULIC,LSC),col="black",lty=2)

text("Phase I",x=3,y=80)
text("Phase II",x=28,y=80)


legend("topright",legend=c("Proposed Method","Naive Method"),lty=c(1,2),
       col=c("red", "blue" ))

