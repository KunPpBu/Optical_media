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
                  x2=x2_stan)
# Load Stan file
stan_code <- "New_Stan.stan"

# Run Stan
fitStan <- stan("New_Stan.stan",data=stan_data, 
                chains = 3, iter = 8500, warmup = 1000, thin = 100, init_r = .1)
print(fitStan, pars=c("sigma"))
plot(fitStan)

pdf(paste(outdir,"trace plot_B.pdf",sep=""))
rstan::traceplot(fitStan, pars = c("B"), inc_warmup = FALSE)
dev.off()

