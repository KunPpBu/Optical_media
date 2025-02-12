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
# t_ijk <- ISO_data$Hours  # Predictors
t_ijk <- rep(seq(100,2500,500),90)
y_ijk<- ISO_data$`Deg. data`
T <- ISO_data$Temp.
RH <- ISO_data$HR
logRH <- log(RH,base=exp(1))
#set the transformed data x1 and x2
x1_stan <- (logRH-log(40))/(log(85)-log(40))
head(x1_stan)
x2_stan <- (11605/(T+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15))
head(x2_stan)

# Set up the initial value of parameters
mu <- rnorm(450,1,1)
sigma <- rnorm(450,7,2)
mat <- matrix(runif(450,1,3),nrow=N, ncol=2)
mu_mat <- matrix(runif(450,1,3),nrow=N,ncol=2)
sigma_mat <- matrix(runif(450,0,2),nrow=2, ncol=2)
A <- rnorm(N,10,3)
B <- rnorm(N,10,3)
delta_H <- rnorm(450,10,1.5)
sigma_for_prior <- matrix(c(1,1,1,1),2,2)
# beta1<- exp(log(A) + B*x1_stan + delta_H*x2_stan)
# for(n in 1:N){
#   for(k in 1:K){
#     beta1[(n-1)*5+k] <- exp(log(A[n]) + B[n]*x1_stan[(n-1)*5+k] + delta_H[n]*x2_stan[(n-1)*5+k]) 
#   }
# }

# beta0 <- mat[,1]
# gamma <- mat[,2]

# Set up the initial value of outcome
for(n in 1:N){
  beta0[n] <- mat[n,1]
  gamma[n] <- mat[n,2]
  for(k in 1:K){
    beta1[(n-1)*5+k] <- exp(log(A[n]) + B[n]*x1_stan[(n-1)*5+k] + delta_H[n]*x2_stan[(n-1)*5+k])
    mu[(n-1)*5+k] <- beta0[n] + beta1[(n-1)*5+k] * (t_ijk[(n-1)*5+k]^gamma[n])
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
stan_code <- "New_Stan.stan"

# cat(stan_code)

# Run Stan
runStan <- stan("New_Stan.stan",data=stan_data, 
                chains = 3, iter = 8500, warmup = 1000, thin = 100, init_r = .1)
print(runStan, pars=c("mat"))


list_of_draws <- extract(runStan)
print(names(list_of_draws))
head(list_of_draws$sigma_mat)
dim(list_of_draws$mat)







# saveRDS(runStan, "runstan.rds")
print(runStan, pars=c("mat","A","B","delta_H","sigma_mat","sigma"))


resStanExt <- rstan::extract(runStan, permuted = TRUE)
# rstan::traceplot(runStan, pars = c("mat","A","B","delta_H","sigma_mat","sigma"), inc_warmup = FALSE)
# rstan::traceplot(runStan, pars = c("A"), inc_warmup = FALSE)
pdf(paste(outdir,"stan trace plot_A.pdf",sep=""))
rstan::traceplot(runStan, pars = c("A"), inc_warmup = FALSE)
dev.off()



