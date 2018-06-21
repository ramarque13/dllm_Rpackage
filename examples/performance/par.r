
  rm(list=ls())
 
  library(dllmsmc); require(mgcv); require(smcUtils)
  
 
  #****************** Particle Parameter Estimation for  *****************************#

    #****************** Stochastic Volatility Models  *****************************#


library(dllmsmc)

sigma.x <- 1
phi <- 0.9
alpha <- -3	 
beta0 <- alpha*(1-phi)
tau.x <- 1/sigma.x^2
tau.inla = tau.x*(1-phi^2)

nT <- 25000
tstart.bugs = nT - 1000
N <- 2000
burnin <- N
tstart = 250  # Number of MCMC steps within PF and Blocking 
nthin <- ifelse(N<500,10,2)
thin.mcmc <- 100
b <- 50
tau0 <- 0.001
b.initial = 2500  # the time when we start the blocking  
N.mcmc = 3000  

S <- 1

D.alpha <- array(NA,c(S,nT,2))
D.phi <- array(NA,c(S,nT,2))
D.tau <- array(NA,c(S,nT,2))

set.seed(106)

#Simulate data
d <- dllm.sim(nT,alpha=alpha,sigma.x=sigma.x,sigma.y=sigma.y,obsmodel="svmodel",phi=phi)

plot.ts(d$x)
plot.ts(d$y)


            #  Unkown static Parameters #

alphafixed <- F
phifixed <- F
taufixed <- F

seed=    sample(1:10^6,1)
                   
res.bugs = dllm.bugs(d$y,obsmodel="svmodel", ,
                     alphaprior=list(mu=0,prec=1,initial=beta0,fixed=alphafixed),
                     phiprior=list(mu=0,prec=1,initial=phi,fixed=phifixed),
                     tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
                     N=N, burnin=burnin,nthin=nthin)                    
    
setwd("folder.to.save.the.outputs")
write.table(data.frame(res.bugs$alpha,res.bugs$phi,res.bugs$tau,
res.bugs$meanT,res.bugs$meanSQ,res.bugs$meanC),paste("res.bugs_nT=",nT,".txt",sep=""))
    
rm(res.bugs)
    
#Ordinary SMC
res.smc <- dllm.smc(d$y,obsmodel="svmodel",initial.method="bugs",
                    alphaprior=list(mu=0,prec=1,initial=beta0,fixed=alphafixed),
                    phiprior=list(mu=0,prec=1,initial=phi,fixed=phifixed),
                    tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
                    N=N,start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),tau0=tau0,seed=seed)

write.table(data.frame(res.smc$alpha[nT,],res.smc$phi[nT,],res.smc$tau[nT,],
res.smc$meanT[nT,],res.smc$meanSQ[nT,],res.smc$meanC[nT,]),paste("res.smc_nT=",nT,".txt",sep=""))
    
    
rm(res.smc)

                    
#res.smc.dyn <- summary(res.smc)

b=50
res.block <- dllm.smc(d$y,obsmodel="svmodel",initial.method="bugs",
                      alphaprior=list(mu=0,prec=1,initial=beta0,fixed=alphafixed),
                      phiprior=list(mu=0,prec=1,initial=phi,fixed=phifixed),
                      tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
                      N=N,b=b,start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),tau0=tau0,b.initial=2500,k=b,
                      seed=seed)

                      
write.table(data.frame(res.block$alpha[nT,],res.block$phi[nT,],res.block$tau[nT,],
res.block$meanT[nT,],res.block$meanSQ[nT,],res.block$meanC[nT,]),paste("res.smc_nT=",nT,".txt",sep=""))
rm(res.block)

res.bugs = read.table("res.bugs_nT=10000.txt"); names(res.bugs) = c("alpha","phi","tau","meanT","meanSQ","meanC")

res.smc = read.table("res.smc_nT=10000.txt"); names(res.smc) = c("alpha","phi","tau","meanT","meanSQ","meanC")
res.block = read.table("res.block_nT=10000.txt"); names(res.block) = c("alpha","phi","tau","meanT","meanSQ","meanC")


par.col = c("black","gray","navyblue")

par(mfrow=c(1,1))
d1 <- density(res.bugs$meanT)
d2 <- density(res.smc$meanT)
d3 <- density(res.block$meanT)
matplot(cbind(d1$x,d2$x,d3$x),cbind(d1$y,d2$y,d3$y),type="l",lty=1,main="",col=par.col,lwd=2.5,xlab="",ylab="",cex.axis=1.25)

par(mfrow=c(1,1))
d1 <- density(res.bugs$meanSQ)
d2 <- density(res.smc$meanSQ)
d3 <- density(res.block$meanSQ)
matplot(cbind(d1$x,d2$x,d3$x),cbind(d1$y,d2$y,d3$y),type="l",lty=1,main="",col=par.col,lwd=2.5,xlab="",ylab="",cex.axis=1.25)

par(mfrow=c(1,1))
d1 <- density(res.bugs$meanC)
d2 <- density(res.smc$meanC)
d3 <- density(res.block$meanC)
matplot(cbind(d1$x,d2$x,d3$x),cbind(d1$y,d2$y,d3$y),type="l",lty=1,main="",col=par.col,lwd=2.5,xlab="",ylab="",cex.axis=1.25)

par(mfrow=c(1,1))
d1 <- density(res.bugs$alpha)
d2 <- density(res.smc$alpha)
d3 <- density(res.block$alpha)
matplot(cbind(d1$x,d2$x,d3$x),cbind(d1$y,d2$y,d3$y),type="l",lty=1,main="",col=par.col,lwd=2.5,xlab="",ylab="",cex.axis=1.25)
abline(v=alpha,lty=3,col="red",lwd=2.5)

par(mfrow=c(1,1))
d1 <- density(res.bugs$phi)
d2 <- density(res.smc$phi)
d3 <- density(res.block$phi)
matplot(cbind(d1$x,d2$x,d3$x),cbind(d1$y,d2$y,d3$y),type="l",lty=1,main="",col=par.col,lwd=2.5,xlab="",ylab="",cex.axis=1.25)
abline(v=phi,lty=3,col="red",lwd=2.5)
#legend("topright",legend=c("Bugs","SMC","BLOCK"),lty=1,col=1:3)
dev.off()

par(mfrow=c(1,1))
d1 <- density(res.bugs$tau)
d2 <- density(res.smc$tau)
d3 <- density(res.block$tau)
matplot(cbind(d1$x,d2$x,d3$x),cbind(d1$y,d2$y,d3$y),type="l",lty=1,main="",col=par.col,lwd=2.5,xlab="",ylab="",cex.axis=1.25)
abline(v=tau.x,lty=3,col="red",lwd=2.5)
