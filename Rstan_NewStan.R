#load libraries
library(ggplot2)
library(StanHeaders)
library(rstan)
library(RColorBrewer)
#where the STAN model is saved
setwd("~/GitHub/Optical_media/")
## Data generation
set.seed(3117)
#set up the model data
N <- 90
K <- 5
t_ijk <- ISO_data$Hours  # Predictors
y_ijk<- ISO_data$`Deg. data`
T <- ISO_data$Temp.
RH <- ISO_data$HR
logRH <- log(RH,base=exp(1))
#set the transformed data x1 and x2
x1_stan <- (logRH-log(40))/(log(85)-log(40))
head(x1_stan)
x2_stan <- (11605/(T+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15))
head(x2_stan)

#set the parameters
beta0 <- 3
beta1 <- 3
sigma <- 1.5
mat <- matrix(runif(450,0,30),nrow=N, ncol=2)
mu_mat <- matrix(runif(450,0,3),nrow=N,ncol=2)
sigma_mat <- matrix(runif(450,0,2),nrow=2, ncol=2)
A <- runif(90,70,80)
B <- runif(90,-3,3)
delta_H <- runif(90,0.5,1.5)
sigma_for_prior <- diag(1,2,2)
for(n in 1:N){
  for(k in 1:K){
    beta1[(n-1)*5+k] <- exp(log(A[n]) + B[n]*x1[(n-1)*5+k] + delta_H[n]*x2[(n-1)*5+k]) 
  }
}

beta0 <- mat[,1]
gamma <- mat[,2]
mu <- beta0 + beta1 * (t_ijk^mat[,2])
# outcome

for(n in 1:N){
  for(k in 1:K){
    y_ijk[(n-1)*5+k] <- mu[(n-1)*5+k] + sigma[(n-1)*5+k]
  }
}

# Set model data
stan_data <- list(N=N,
                  K=K,
                  t_ijk=t_ijk,
                  y_ijk=y_ijk,
                  RH=RH,
                  T=T,
                  x1=x1_stan,
                  x2=x2_stan)

# Load Stan file
fileName <- "New_Stan.stan"
stan_code <- readChar(fileName,file.info(fileName)$size)
cat(stan_code)

# Run Stan
runStan <- stan(model_code=stan_code,data=stan_data,
                chains = 3, iter = 3000,warmup = 500, thin = 10)






