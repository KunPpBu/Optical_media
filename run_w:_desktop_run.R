load("~/GitHub/Optical_media/stan_sampler 2.RData")
ss_desktop_run<- get(load("~/GitHub/Optical_media/stan_sampler 2.RData"))
dim(ss_desktop_run)
library("mvtnorm")
library(MASS)
N <- 29850
mu <- rep(0,2)

sigma <- list()
bvn <- list()

for(n in 1:N){
  sigma[[n]] <- matrix(c(ss_desktop_run[n,7],ss_desktop_run[n,8],ss_desktop_run[n,8],ss_desktop_run[n,9]),nrow=2,ncol=2) 
  bvn[[n]] <- mvrnorm(1, mu = mu, Sigma = sigma[[n]] )
}
bvn<-matrix(unlist(bvn),ncol=2,byrow=T)
#set a time 
t <- seq(0,109999,by=1000)
#input the normal use condition x1 is RH, x2 is Temp
x1 <- rep((log(50)-log(40))/(log(85)-log(40)),110)
head(x1)
x2 <- rep((11605/(25+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15)),110)
head(x2)

#function for degrad
re_func <- function(t, par, bvn,x1,x2){
  beta0 <- NULL
  beta1 <- NULL
  gamma <- NULL
  deg <- matrix(rep(0,N*110),N,110)
  for(i in 1:N){
    for(j in 1:110)
    {
      beta0[i] <- par[i,4]+bvn[i,1]
      beta1[i] <- exp(log(par[i,3])+par[i,1]*x1[j]+par[i,2]*x2[j])
      gamma[i] <- par[i,5]+bvn[i,2]
      deg[i,j] <- beta0[i]+beta1[i]*(t[j]^(gamma[i]))
    }
  }
  return(deg)
}
deg<-re_func(t=t,par=ss_desktop_run,bvn=bvn,x1=x1,x2=x2)
deg_prob<-(deg>log(280))
deg_new<-apply(deg_prob,2,mean)
# plot(1-deg_new)
library(zoom)
pdf(paste(outdir,"re_func_with_ss_desktop.pdf",sep=""))
plot(1-deg_new)
dev.off()
