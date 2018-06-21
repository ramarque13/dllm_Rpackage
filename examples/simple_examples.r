
  # Simple examples
  
  sigma.x <- 1
  sigma.y <- 0.5
  phi <- 0.5
  alpha <- 0
  beta0 <- alpha*(1-phi)
  tau.x <- 1/sigma.x^2
  tau.y <- 1/sigma.y^2
  nT <- 500

data.gaus <- dllm.sim(nT,alpha=alpha,sigma.x=sigma.x,sigma.y=sigma.y, 
obsmodel="lingauss", phi=phi)
plot.ts(data.gaus\$y,xlab="",ylab="")

data.sv <- dllm.sim(nT,alpha=alpha,sigma.x=sigma.x,sigma.y=sigma.y, 
obsmodel="svmodel", phi=phi)
plot.ts(data.sv\$y,xlab="",ylab="")

N <- 1000
tstart <- 25
b <- 50 # block size
alphafixed <- phifixed <- taufixed <- TRUE
seed <- sample(1:9999,1)
 #Ordinary SMC
  res.smc <- dllm.smc(data.sv$y,obsmodel="svmodel",
                      alphaprior=list(mu=0,prec=1e-04,initial=beta0,fixed=alphafixed),
                      phiprior=list(mu=0,prec=1e-04,initial=phi,fixed=phifixed),
                      tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
                      N=N,start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),tau0=tau0,seed=seed)
  
 res.smc.dyn <- summary(res.smc) 
 plot(res.smc.dyn) 
 
 #Blocking SMC
 b=50  #blocklength
 res.block <- dllm.smc(data.sv$y,obsmodel="svmodel",
                       alphaprior=list(mu=0,prec=1e-04,initial=beta0,fixed=alphafixed),
                       phiprior=list(mu=0,prec=1e-04,initial=phi,fixed=phifixed),
                       tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
                       N=N,b=b,start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),k=b,
                       seed=seed)
                       
 res.block.dyn <- summary(res.block)
 plot(res.block.dyn)

b=100  #blocklength
res.block2 <- dllm.smc(data.sv$y,obsmodel="svmodel",
                       alphaprior=list(mu=0,prec=1e-04,initial=beta0,fixed=alphafixed),
                       phiprior=list(mu=0,prec=1e-04,initial=phi,fixed=phifixed),
                       tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
                       N=N,b=b,start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),k=b,
                       seed=seed)
                       
res.block.dyn2 <- summary(res.block2)
plot(res.block.dyn,res.block.dyn,res.block.dyn2)

alphafixed <- phifixed <- taufixed <- FALSE

res.block <- dllm.smc(data.sv$y,obsmodel="svmodel",
                       alphaprior=list(mu=0,prec=1e-04,initial=beta0,fixed=alphafixed),
                       phiprior=list(mu=0,prec=1e-04,initial=phi,fixed=phifixed),
                       tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
                       N=N,b=b,start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),k=b,
                       seed=seed)
                       
res.block.dyn <- summary(res.block)
  
t = 300
par.block <- data.frame(res.block.dyn\$alpha[t,],res.block.dyn\$phi[t,],res.block.dyn\$tau[t,])
colnames(par.block) <- c("alpha","phi","tau") 
par.block  
