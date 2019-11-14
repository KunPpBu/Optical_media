#load libraries
library(ggplot2)
library(StanHeaders)
library(rstan)
library(RColorBrewer)
#where the STAN model is saved
setwd("~/GitHub/Optical_media/")

## Data generation
set.seed(3117)
N <- 90
K <- 5
t_ijk <- rep(seq(100,2500,500),90)
y_ijk<- ISO_data$`Deg. data`
T <- ISO_data$Temp.
RH <- ISO_data$HR
logRH <- log(RH,base=exp(1))
mean_mat <- as.vector(c(0,0),mode = "any")

#set the transformed data x1 and x2
x1_stan <- (logRH-log(40))/(log(85)-log(40))
head(x1_stan)
x2_stan <- (11605/(T+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15))
head(x2_stan)

# Set model data
stan_data <- list(N=N,
                  K=K,
                  t_ijk=t_ijk,
                  y_ijk=y_ijk,
                  RH=RH,
                  T=T,
                  x1=x1_stan,
                  x2=x2_stan,
                  mean_mat=mean_mat)

#initial values
sigma_for_prior <- matrix(c(1,1,1,1),2,2)
sigma <- rnorm(450,7,2)
sigma_mat <- matrix(c(1,0,0,1),nrow=2, ncol=2)
A <- rnorm(N,100,3)
B <- rnorm(N,10,3)
delta_H <- rnorm(450,10,1.5)

# Set model code
# Load Stan file
stan_code <- "Stan_test3.stan"

# Run Stan
fitStan <- stan("Stan_test3.stan",data=stan_data, 
                chains = 3, iter = 3500, warmup = 1000, thin = 100, init_r = .1)
print(fitStan, pars=c("A"))
print(fitStan)
plot(fitStan,pars=c("A","B","delta_H"))

# Retrieve the posterior draws 
matrix_of_draws <- as.matrix(fitStan)
print(colnames(matrix_of_draws))
sampler_draws <- as.matrix(fitStan, pars = c("mu_mat[1]", "mu_mat[2]","sigma","A","B","delta_H","sigma_mat[1,1]","sigma_mat[1,2]","sigma_mat[2,2]"))
head(sampler_draws)
mean(sampler_draws[,1])
# Plot parameters trace #individually
pdf(paste(outdir,"trace plot_ABH1.pdf",sep=""))
rstan::traceplot(fitStan, pars = c("A","B","delta_H"), inc_warmup = FALSE)
dev.off()
# Plot all parameters trace in one pdf file
pdf(paste(outdir,"trace plot_newstan1.pdf",sep=""))
rstan::traceplot(fitStan, pars = c("mu_mat[1]","mu_mat[2]","delta_H","A","B","sigma","sigma_mat[1,1]","sigma_mat[2,2]"), inc_warmup = FALSE)
dev.off()

