#inputs for the function

#mcmc parameters simulation result
pmat<- get(load("~/Documents/Kun Bu/simulation/coda_new.rdata")) 
dim(pmat)

test <- get(load("~/Documents/Kun Bu/simulation/coda_new.rdata"))
dim(test)

#generate the random effect 
library("mvtnorm")
library(MASS)
N <- 3000
mu <- rep(0,2)

sigma <- list()
bvn <- list()

for(n in 1:3000){
  sigma[[n]] <- matrix(c(pmat[n,7],pmat[n,8],pmat[n,8],pmat[n,10]),nrow=2,ncol=2) 
  bvn[[n]] <- mvrnorm(1, mu = mu, Sigma = sigma[[n]] )
}
bvn<-matrix(unlist(bvn),ncol=2,byrow=T)
#set a time 
t <- seq(0,10599999,by=2000)
#input the normal use condition x1 is RH, x2 is Temp
x1 <- rep((log(50)-log(40))/(log(85)-log(40)),53000)
head(x1)
x2 <- rep((11605/(25+273.15)-11605/(15+273.15))/(11605/(85+273.15)-11605/(15+273.15)),5300)
head(x2)

#function for degrad
re_func <- function(t, par, bvn,x1,x2){
  beta0 <- NULL
  beta1 <- NULL
  gamma <- NULL
  deg <- matrix(rep(0,3000*53000),3000,5300)
  for(i in 1:3000){
    for(j in 1:5300)
    {
      beta0[i] <- par[i,4]+bvn[i,1]
      beta1[i] <- exp(par[i,3]+par[i,1]*x1[i]+par[i,2]*x2[i])
      gamma[i] <- par[i,5]+bvn[i,2]
      deg[i,j] <- beta0[i]+beta1[i]*(t[j]^(gamma[i]))
    }
  }
  return(deg)
}
deg<-re_func(t=t,par=test,bvn=bvn2_,x1=x1,x2=x2)
deg_prob<-(deg>log(280))
deg_new<-apply(deg_prob,2,mean)
# plot(1-deg_new)
library(zoom)
pdf(paste(outdir,"re_func.pdf",sep=""))
plot(1-deg_new)
dev.off()






ggplot() + 
  geom_point(mapping = aes(x = t[1:800], y = deg[1:800]))
ggplot()+
  geom_smooth(mapping=aes(x=t[1:800],y=1-deg_new[1:800]))  

