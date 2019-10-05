ISO_data <- read_excel("~/Google Drive/School/USF/Spring 2019/Master's Thesis/Optical Media/ISO data.xlsx")
RH <- ISO_data$HR
logRH <- log(RH,base=exp(1))
x1 <- (logRH-log(40))/(log(85)-log(40))
T <- ISO_data$Temp.
x2 <- (11605/(T+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15))
t <- ISO_data$Hours
y <- ISO_data$`Deg. data`
N <- 90
M <- 5


model_code=cat('model
               {
               for (i in 1:N){
               mat[i,1:2] ~ dmnorm(mu_mat[,],prec.sigma[,])
               for (k in 1:M){
               y[(i-1)*5+k] ~ dlnorm (mu[(i-1)*5+k], tau ) 
               mu[(i-1)*5+k] <- mu_mat[,1] + bvn[,1] + (beta1[(i-1)*5+k]*pow(t[(i-1)*5+k],mu_mat[,2]+bvn[,2]))
               beta1[(i-1)*5+k] <- exp(logA + B*x1[(i-1)*5+k] + delta_H*x2[(i-1)*5+k])
               }
               }
               
               #priors
               mu_mat[1,1] ~ dnorm(0,1e-3)
               mu_mat[1,2] ~ dnorm(0,1e-3)
               bvn[1,1] ~ dnorm(0,1e-3)
               bvn[1,2] ~ dnorm(0,1e-3)
               prec.sigma[1:2,1:2] ~ dwish(R[1:2,1:2],2)
               sigma1[1:2,1:2] <- inverse(prec.sigma[,])
               rho12 <- sigma1[1,2]/sqrt(sigma1[1,1]*sigma1[2,2])
               R[1,1] <- 1
               R[1,2] <- 0
               R[2,1] <- 0
               R[2,2] <- 1
               logA ~ dnorm(0,1e-4)
               B ~ dnorm(0, 1e-2)
               delta_H ~ dnorm(0,1e-4)
               tau ~ dgamma(1e-3,1e-3)           
               sigma <- 1/sqrt(tau)
               
               }', file="lognormal.jags"
)

model_parameters =  c( "logA", "B", "delta_H","mu_mat","bvn", "sigma","sigma1")
model.data <- list(x1=x1,x2=x2,t=t,y=y,N=N,M=M)

#install.packages('rjags')
library('rjags')
model_run <- jags.model(file="lognormal.jags",
                        data = model.data,
                        
                        n.chains = 3,
                        n.adapt = 50000)

samples.jags <- jags.samples(model_run,
                             model_parameters,
                             n.iter=20000, 
                             thin=6500,
                             n.adapt=500000)

samples.coda <- coda.samples(model_run,
                             model_parameters,
                             n.iter=200000,
                             thin=6500,
                             n.adapt=500000
)

summary(samples.coda)
pdf(paste(outdir,"trace plot_bvn.pdf",sep=""))
plot(samples.coda[[1]][,1:12]) 
dev.off()  
# save the mcmc results
pest1=summary(samples.coda)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
pmat1=as.matrix(samples.coda)
outdir=""
save(pmat1,file=paste(outdir,"coda_with_random_effect",".rdata",sep=""))
#ACF plot
pdf(paste(outdir,"ACF_bvn.pdf",sep=""))
par(mfrow=c(2,2))
acf(pmat,xlab="Lag",ylab="Correlation") # to improve acf plot, thin=1000, used to be 100
dev.off()