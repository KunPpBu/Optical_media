#load libraries

library(ggplot2)
library(StanHeaders)
library(rstan)
library(RColorBrewer)
#where the STAN model is saved
setwd("~/Documents/Kun_Bu/simulation/")
set.seed(1000000)
#set up the model data
N <- 90
K <- 5
t <- ISO_data$Hours
y <- ISO_data$`Deg. data`
T <- ISO_data$Temp.
RH <- ISO_data$HR
logRH <- log(RH,base=exp(1))
id<-rep(1:K,each=(N-1)*5+K) #index of plant species

#standard deviation of the group-level coefficient
sigma_for_prior<- matrix(c(1,2,2,1),nrow=N,ncol=N)
sigma_mat1 <- matrix(c(1,1,1,1), nrow=N, ncol=N)
sigma_mat <-as.numeric(sigma_mat1)

#set up for the initial value of parameters
tau<-c(0.3,2,1)
sigma<-0
logA <- 0
B <- 0
delta_H <-0
beta1<-0 # beta1 is the transformed parameters
mat<-matrix(c(1,1),nrow=N,ncol=N)

mat<-1
mu <- 1
mu_mat1 <- matrix(c(0,0),nrow = N,ncol = N)
mu_mat <- as.numeric(mu_mat1) #mu_mat[,1] is mu0, mu_mat[,2] is gamma0


#set the transformed data x1 and x2
x1_stan <- (logRH-log(40))/(log(85)-log(40))
head(x1_stan)
x2_stan <- (11605/(T+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15))
head(x2_stan)

#group-level regression coefficients
# beta1[(N-1)*5+K]<-exp(logA + B*x1[(N-1)*5+K] + delta_H*x2[(N-1)*5+K]);
# mu[(N-1)*5+K] = mat[N,1] + (beta1[(N-1)*5+K]*((t[(N-1)*5+K])%^%(mat[N,2])));
#
# #the model matrix
# X<-model.matrix(~x+y,data=data.frame(x=runif(N,-2,2),y=runif(N,-2,2)))

# set up the model
for(n in 1:N){
  #simulate response data
  mat[n,1:2] ~ dmnorm(mu_mat[,],sigma[,])
  for(k in 1:K){
    y[(n-1)*5+k] ~ rnorm(mu[(n-1)*5+k], sigma)
    mu[(n-1)*5+k] <- mat[n,1] + (beta1[(n-1)*5+k]*(t[(n-1)*5+k]^mat[n,2]))
    beta1[(n-1)*5+k] <- exp(logA + B*x1[(n-1)*5+k] + delta_H*x2[(n-1)*5+k])
  }
  return(y)}

#run the model--initialization
data_stan <- list(N=N,K=K,x1=x1_stan,x2=x2_stan,y=y,t=t,RH=RH,wishart_df=100)
m_hier<-stan(file="Rstan optical media.stan",data=data_stan,cores = 4, chains = 1, warmup = 250, iter = 500)

