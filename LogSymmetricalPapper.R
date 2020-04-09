############### Author: LUCAS SALES ############
#Institution: Universidade Federal do Rio Grande do Norte
#e-mail: lucasofsales@gmail.com  or ldo-sales@hotmail.com
#Log-symmetrical Class simulation #######


### This code was developed for distribution which belongs to log-symmetrical
#clas. More especifically, the code was illustrated for the log-normal distribution.
#However, to use other distributions of this class, change the rlnorm 
#to the desired distribution.


library(psych)
library(stats)
library(ICSNP)
library(xtable)
library(LaplacesDemon) 

#Proposed estimator I - Balakrishnan et al 2017
PEI=function(x){
  sqrt(mean(x)*harmonic.mean(x)) 
}


#n is the sample size
#med is the median of the distribution
#delta is the shif in the median
#q1 and q2 is the control limits
#d is the bias correct factor of the estimator.

#Generate PEI samples
amostraPEI=function(n,med,delta){
  x=NULL
  for(i in 1:10000){
    # print(i)
    x[i]=PEI(rlnorm(n,meanlog=log(med*(1+delta)))) #log-normal distribution
  }
  return(x)
}
#Estimated ARL of PEI
arlPEI=function(n,med,delta,q1,q2){
  x1=amostraPEI(n,med,delta)
  z1=(x1-med)/sd(x0)
  fda1=ecdf(z1)
  (arl1=1/(fda1(q1)+1-fda1(q2)))   
  (arl1f=1/(fda1(-3)+1-fda1(3)))
  return(c(arl1,arl1f))
}

#Generate PEIII (hodgers-lehmann estimator) samples
amostraHL=function(n,med,d,delta){
  x=NULL
  for(i in 1:10000){
    # print(i)
    x[i]=hl.loc(rlnorm(n,meanlog=log(med*(1+delta))))/(1+d) #log-normal distribution
  }
  return(x)
}
#Estimated ARL of PEIII
arlHL=function(n,med,d,delta,q1,q2){
  x1=amostraHL(n,med,d,delta)
  z1=(x1-med)/sd(x0)
  fda1=ecdf(z1)
  (arl1=1/(fda1(q1)+1-fda1(q2)))   
  (arl1f=1/(fda1(-3.09)+1-fda1(3.09)))
  return(c(arl1,arl1f))
}

#Example of how to obtain the control limits of a log-normal distribution
#using the PEI. 

med=2
n=10000
alpha=0.00198

x0=amostraPEI(n,med,0)
z0=(x0-med)/sd(x0)
(q1=quantile(z0,probs=c(alpha/2)))
(q2=quantile(z0,probs=c(1-alpha/2)))
fda0=ecdf(z0)
