#load libraries

library(ggplot2)
library(rstan)
library(StanHeaders)
library(RColorBrewer)
#where the STAN model is saved
setwd("~/Documents/Kun Bu/simulation/")
set.seed(1000000)
N <- 90
K <- 5
R <- matrix(c(1,2,2,1),nrow=N,ncol=N)
id<-rep(1:K,each=(N-1)*5+K) #index of plant species
#standard deviation of the group-level coefficient
tau<-c(0.3,2,1)
sigma<-0
logA <- 0
B <- 0
delta_H <-0
beta1<-0
mat<-matrix(NA,nrow=N,ncol=2)
mu <- 0
mu_mat <- matrix(NA,nrow = N,ncol = 2)

x1_stan <- rep((log(50)-log(40))/(log(85)-log(40)),450)
head(x1)
x2_stan <- rep((11605/(25+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15)),450)
head(x2)

#group-level regression coefficients
beta1[(N-1)*5+K]<-exp(logA + B*x1[(N-1)*5+K] + delta_H*x2[(N-1)*5+K]);
# mu[(N-1)*5+K] = mat[N,1] + (beta1[(N-1)*5+K]*((t[(N-1)*5+K])%^%(mat[N,2])));
#the model matrix
X<-model.matrix(~x+y,data=data.frame(x=runif(N,-2,2),y=runif(N,-2,2)))
y <- ISO_data$`Deg. data`
length(y)

for(n in 1:N){
  #simulate response data
  mat[n,1:2] ~ dmnorm(mu_mat[,],sigma[,])
  for(k in 1:K){
  y[(n-1)*5+k] ~ rnorm(mu[(n-1)*5+k], sigma)
  
}}
#run the model--initialization
m_hier<-stan(file="Rstan optical media.stan",data=list(N=N,K=K,x1=x1_stan,x2=x2_stan,y=y,t=t))

