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
mu_mat <- matrix(runif(450,1,2),nrow=N,ncol=2)
sigma_for_prior <- matrix(c(1,1,1,1),2,2)
sigma <- rnorm(450,0,1)
sigma_mat <- matrix(c(1,0.5,0.5,1),nrow=2, ncol=2)
A <- rnorm(N,0,2)
B <- rnorm(N,0,2)
delta_H <- rnorm(450,10,1.5)
# epsilon <- matrix(runif(450,1,2),nrow=N,ncol=2)

# Set model code
# Load Stan file
stan_code <- "Stan_test3.stan"

# Run Stan
fitStan <- stan("Stan_test3.stan",data=stan_data, 
                chains = 3, iter = 5500, warmup = 1000, thin = 100, init_r = .1)
print(fitStan, pars=c("A"))
plot(fitStan,pars=c("A","B","delta_H","sigma","sigma_mat","mu_mat"))

# Retrieve the posterior draws 
matrix_of_draws <- as.matrix(fitStan)
print(colnames(matrix_of_draws))
sampler_draws <- as.matrix(fitStan, pars = c("B", "delta_H","A","mu_mat[1]","mu_mat[2]","sigma","sigma_mat[1,1]","sigma_mat[1,2]","sigma_mat[2,2]"))
head(sampler_draws)
dim(sampler_draws)
mean(sampler_draws[,8])
# Save Stan Output
save(sampler_draws, file = "stan_sampler.RData")

# mcmc_pairs
model_cp<-stan_model("Stan_test3.stan")
fit_cp <- sampling(model_cp, data = stan_data, seed = 803214053, control = list(adapt_delta = 0.9))
posterior_cp <- as.array(fit_cp)
lp_cp <- log_posterior(fit_cp)
np_cp <- nuts_params(fit_cp)
mcmc_pairs(posterior_cp, np = np_cp, pars = c("B", "delta_H","A","mu_mat[1]","mu_mat[2]","sigma"),
           off_diag_args = list(size = 0.75))
color_scheme_set("red")
mcmc_nuts_divergence(np_cp, lp_cp)
color_scheme_set("red")
mcmc_nuts_energy(np_cp)
ratios_cp <- neff_ratio(fit_cp)
mcmc_neff(ratios_cp, size = 2)


# Plot parameters trace #individually
pdf(paste(outdir,"trace plot_ABH1.pdf",sep=""))
rstan::traceplot(fitStan, pars = c("A","B","delta_H"), inc_warmup = FALSE)
dev.off()
# Plot all parameters trace in one pdf file
pdf(paste(outdir,"trace plot_newstan_test4.pdf",sep=""))
rstan::traceplot(fitStan, pars = c("mu_mat[1]","mu_mat[2]","delta_H","A","B","sigma","sigma_mat[1,1]","sigma_mat[1,2]", "sigma_mat[2,2]"), inc_warmup = FALSE)
dev.off()

# ACF plots
install.packages("bayesplot")
library("bayesplot")
pdf(paste(outdir,"ACF_test4.pdf",sep=""))
posterior_cp <- as.array(fitStan)
mcmc_acf(posterior_cp, pars = c("A","B","delta_H","sigma","mu_mat[1]","mu_mat[2]","sigma_mat[1,1]","sigma_mat[1,2]", "sigma_mat[2,2]"), lags = 20)
dev.off()

