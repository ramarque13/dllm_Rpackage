


  rm(list=ls())
 
  library(dllmsmc); require(mgcv); require(smcUtils)
  
 
  #****************** Simulating the path functionals *****************************#

  prin <- FALSE
  #Note: This script is time-consuming!!

  MC.runs = 50
  sigma.x <- 1
  sigma.y <- 1
  phi <- 0.9
  alpha <- 0
  beta0 <- alpha*(1-phi)
  tau.x <- 1/sigma.x^2
  tau.y <- 1/sigma.y^2
  nT <- 30000
  N <- 1000
  burnin <- N
  tstart <- 10
  nthin <- 1
  thin.mcmc <- 50
  b <- 50*1
  tau0 <- 0.001
  S <- 1

  D.alpha <- array(NA,c(S,nT,2))
  D.phi <- array(NA,c(S,nT,2))
  D.tau <- array(NA,c(S,nT,2))

  set.seed(11) #83
  seed = 1234

  #dataset
  d <- dllm.sim(nT,alpha=alpha,sigma.x=sigma.x,sigma.y=sigma.y,obsmodel="lingauss",phi=phi)
  y <- d$y; rm(d)
  # Save dir
  save.path.data = "folder.to.save.the.simulated.dataset"

   # fixed parameters 
  alphafixed <- T;  phifixed <- T;  taufixed <- T
  
  
################################################################
##                                                             ##
##                    Ordinary SMC                             ##
#################################################################
 
 
  #Ordinary SMC
  for(j in 1:MC.runs)
    {    
  res.smc <- dllm.smc(y,obsmodel="lingauss",
		alphaprior=list(mu=0,prec=1e-04,initial=beta0,fixed=alphafixed),
		phiprior=list(mu=0,prec=1e-04,initial=phi,fixed=phifixed),
		tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
		N=N,start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),tau0=tau0,seed=j)

  name.file= paste("data1_smc_N=",N,"_MC",j,".rds",sep="")
  if(N==1000)
    save.dir.res =  paste(save.path.data,"SMC/N=1000/",sep="")
 
  setwd(save.dir.res)		
  saveRDS(res.smc,name.file)
  
  rm(res.smc)
    }
  
  
################################################################
##                                                             ##
##                    Blocking                                 ##
#################################################################


  
  b=c(50,100,200) # blocks
  
  for(jj in 1:length(b))
    {
    if(N==1000)
        save.dir.res =  paste(save.path.data,"Block/N=1000/",sep="")
    for(j in 1:MC.runs)  # salving as seed = j
        {
    res.block <- dllm.smc(y,obsmodel="lingauss",
               alphaprior=list(mu=0,prec=1e-04,initial=beta0,fixed=alphafixed),
               phiprior=list(mu=0,prec=1e-04,initial=phi,fixed=phifixed),
               tauprior=list(shape=1,rate=1,initial=tau.x,fixed=taufixed),
               N=N,b=b[jj],start.mcmc=list(m=tstart,nthin=nthin,burnin=burnin),tau0=tau0,b.initial=1,k=b[jj],
               seed=j)
  name.file= paste("data1_block_N=",N,"_b=",b[jj],"_MC",j,".rds",sep="")
  
  setwd(save.dir.res)		
  saveRDS(res.block,name.file)
  rm(res.block)
    }
        }
        

    
#################################################################
##                                                              ##
##                    MCMC - kalman smoothing                   ##
#################################################################
  kalman = TRUE
  if(kalman)
	{
  seed=1234
  save.dir.res =  paste(save.path.data,"MCMC/",sep="")
  res.mcmc <- dllm.dyn.mcmc(y,start=tstart,
		alphaprior=list(mu=0,prec=0.0001,initial=beta0,fixed=alphafixed),
		phiprior=list(mu=0,prec=0.0001,initial=phi,fixed=phifixed),
		tauprior=list(shape=1,rate=0.1,initial=tau.x,fixed=taufixed),
		N=2000,burnin=burnin,nthin=nthin,thin=thin.mcmc,tau0=tau0,seed=seed)
		
  name.file= paste("data1_MCMC.rds",sep="")
  
  setwd(save.dir.res)		
  saveRDS(res.mcmc,name.file)
  rm(res.mcmc)

  rm(list=ls())
  }
  
 #****************** Time plots of the path functionals *****************************#


save.path.data     = "folder.to.save.the.simulated.dataset"
save.dir.res.smc   =  paste(save.path.data,"SMC/N=",N,"/",sep="")
save.dir.res.block =  paste(save.path.data,"Block/N=",N,"/",sep="")
save.dir.res.mcmc  =  paste(save.path.data,"MCMC/",sep="")
 
 
 
 #meanC = sum x_t x_t-1

 # Ord SMC 

setwd(save.dir.res.smc)
res.smc1      = readRDS(paste0('data1_smc_N=',N,'_MC',1,'.rds'))
res.smc.dyn1  = summary(res.smc1)[1:4]; rm(res.smc1) #"X" "meanT"  "meanSQ" "meanC" 

res.smc2      = readRDS(paste0('data1_smc_N=',N,'_MC',2,'.rds'))
res.smc.dyn2  = summary(res.smc2)[1:4]; rm(res.smc2)

res.smc3      = readRDS(paste0('data1_smc_N=',N,'_MC',3,'.rds'))
res.smc.dyn3  = summary(res.smc3)[1:4]; rm(res.smc3)

res.smc4      = readRDS(paste0('data1_smc_N=',N,'_MC',4,'.rds'))
res.smc.dyn4  = summary(res.smc4)[1:4]; rm(res.smc4)

setwd(save.dir.res.mcmc) 
res.mcmc     = readRDS('data1_MCMC.rds')
res.mcmc.dyn = summary(res.mcmc)[1:4]; rm(res.mcmc)

meanSQ.lim=NULL
meanC.lim=NULL

n <- nrow(res.smc.dyn1$X)
ind1 <- c(1:n)[!is.na(res.smc.dyn1$X[,1])]
ind  <- c(1:n)[!is.na(res.mcmc.dyn$X[,1])]
ind3 <- c(1:n)[!is.na(res.smc.dyn3$X[,1])]
ind4 <- c(1:n)[!is.na(res.smc.dyn4$X[,1])]
ind2 <- c(1:n)[!is.na(res.smc.dyn2$X[,1])]

par.col = c("lightslateblue","gray","black","navyblue","deepskyblue") 

	# posterior mean 
meanC.lim = c(4.5,4.8)
pdf(paste("g1_smc_N=",N,"meanC_post_mean",sep=""))
par(cex.axis=1.6,cex.lab=1.5)
matplot(ind1,res.smc.dyn1$meanC[ind1,1],type="l",lty=1,col=par.col[2],lwd=1.5,xlab="Time",ylab="",ylim=meanC.lim)
matlines(ind,res.mcmc.dyn$meanC[ind,1],type="l",lty=2,lwd=2,col=par.col[3])
matlines(ind3,res.smc.dyn2$meanC[ind3,1],type="l",lty=3,lwd=1.5,col=par.col[2])
matlines(ind4,res.smc.dyn4$meanC[ind4,1],type="l",lty=4,lwd=1.5,col=par.col[2])
matlines(ind2,res.smc.dyn3$meanC[ind2,1],type="l",lty=5,lwd=1.5,col=par.col[2])
dev.off()

	# 2.5% CI
meanC.lim = c(4.5,4.8)
pdf(paste("g1_smc_N=",N,"meanC_post_Q1",sep=""))
par(cex.axis=1.6,cex.lab=1.5)
matplot(ind1,res.smc.dyn1$meanC[ind1,2],type="l",lty=1,col=par.col[2],lwd=1.5,xlab="Time",ylab="",ylim=meanC.lim)
matlines(ind,res.mcmc.dyn$meanC[ind,2],type="l",lty=2,lwd=2,col=par.col[3])
matlines(ind3,res.smc.dyn2$meanC[ind3,2],type="l",lty=3,lwd=1.5,col=par.col[2])
matlines(ind4,res.smc.dyn4$meanC[ind4,2],type="l",lty=4,lwd=1.5,col=par.col[2])
matlines(ind2,res.smc.dyn3$meanC[ind2,2],type="l",lty=5,lwd=1.5,col=par.col[2])
dev.off()

	# 97.5% CI

meanC.lim = c(4.5,4.8)
pdf(paste("g1_smc_N=",N,"meanC_post_Q2",sep=""))
par(cex.axis=1.6,cex.lab=1.5)
matplot(ind1,res.smc.dyn1$meanC[ind1,3],type="l",lty=1,col=par.col[2],lwd=1.5,xlab="Time",ylab="",ylim=meanC.lim)
matlines(ind,res.mcmc.dyn$meanC[ind,3],type="l",lty=2,lwd=2,col=par.col[3])
matlines(ind3,res.smc.dyn2$meanC[ind3,3],type="l",lty=3,lwd=1.5,col=par.col[2])
matlines(ind4,res.smc.dyn4$meanC[ind4,3],type="l",lty=4,lwd=1.5,col=par.col[2])
matlines(ind2,res.smc.dyn3$meanC[ind2,3],type="l",lty=5,lwd=1.5,col=par.col[2])
dev.off()

rm(ind,ind1,ind2,ind3,ind4)
rm(res.mcmc.dyn,res.smc.dyn2,res.smc.dyn3,res.smc.dyn4,res.smc.dyn1)

 # blocking

setwd(save.dir.res.block)
res.block1      = readRDS(paste0('data1_block_N=',N,'_b=',b,'_MC',1,'.rds'))
res.block.dyn1  = summary(res.block1)[1:4]; rm(res.block1) #"X"      "meanT"  "meanSQ" "meanC" 

res.block2      = readRDS(paste0('data1_block_N=',N,'_b=',b,'_MC',2,'.rds'))
res.block.dyn2  = summary(res.block2)[1:4]; rm(res.block2)

res.block3      = readRDS(paste0('data1_block_N=',N,'_b=',b,'_MC',3,'.rds'))
res.block.dyn3  = summary(res.block3)[1:4]; rm(res.block3)

res.block4      = readRDS(paste0('data1_block_N=',N,'_b=',b,'_MC',4,'.rds'))
res.block.dyn4  = summary(res.block4)[1:4]; rm(res.block4)

setwd(save.dir.res.mcmc) 
res.mcmc     = readRDS('data1_MCMC.rds')
res.mcmc.dyn = summary(res.mcmc)[1:4]; rm(res.mcmc)

meanSQ.lim=NULL
meanC.lim=NULL

n <- nrow(res.block.dyn1$X)
ind1 <- c(1:n)[!is.na(res.block.dyn1$X[,1])]
ind  <- c(1:n)[!is.na(res.mcmc.dyn$X[,1])]
ind3 <- c(1:n)[!is.na(res.block.dyn3$X[,1])]
ind4 <- c(1:n)[!is.na(res.block.dyn4$X[,1])]
ind2 <- c(1:n)[!is.na(res.block.dyn2$X[,1])]

 #meanC = sum x_t x_t-1
	  
	# posterior mean 
meanC.lim = c(4.5,4.8)
pdf(paste("g1_block_b=",b,"_N=",N,"meanC_post_mean",sep=""))
par(cex.axis=1.6,cex.lab=1.5)
matplot(ind1,res.block.dyn1$meanC[ind1,1],type="l",lty=1,col=par.col[2],lwd=1.5,xlab="Time",ylab="",ylim=meanC.lim)
matlines(ind,res.mcmc.dyn$meanC[ind,1],type="l",lty=2,lwd=2,col=par.col[3])
matlines(ind3,res.block.dyn2$meanC[ind3,1],type="l",lty=3,lwd=1.5,col=par.col[2])
matlines(ind4,res.block.dyn4$meanC[ind4,1],type="l",lty=4,lwd=1.5,col=par.col[2])
matlines(ind2,res.block.dyn3$meanC[ind2,1],type="l",lty=5,lwd=1.5,col=par.col[2])
dev.off()
	
	# 2.5% CI
meanC.lim = c(4.5,4.8)
pdf(paste("g1_block_b=",b,"_N=",N,"meanC_post_Q1",sep=""))
par(cex.axis=1.6,cex.lab=1.5)
matplot(ind1,res.block.dyn1$meanC[ind1,2],type="l",lty=1,col=par.col[2],lwd=1.5,xlab="Time",ylab="",ylim=meanC.lim)
matlines(ind,res.mcmc.dyn$meanC[ind,2],type="l",lty=2,lwd=2,col=par.col[3])
matlines(ind3,res.block.dyn2$meanC[ind3,2],type="l",lty=3,lwd=1.5,col=par.col[2])
matlines(ind4,res.block.dyn4$meanC[ind4,2],type="l",lty=4,lwd=1.5,col=par.col[2])
matlines(ind2,res.block.dyn3$meanC[ind2,2],type="l",lty=5,lwd=1.5,col=par.col[2])
dev.off()

	# 97.5% CI

meanC.lim = c(4.5,4.8)
pdf(paste("g1_block_b=",b,"_N=",N,"meanC_post_Q2",sep=""))
par(cex.axis=1.6,cex.lab=1.5)
matplot(ind1,res.block.dyn1$meanC[ind1,3],type="l",lty=1,col=par.col[2],lwd=1.5,xlab="Time",ylab="",ylim=meanC.lim)
matlines(ind,res.mcmc.dyn$meanC[ind,3],type="l",lty=2,lwd=3,col=par.col[3])
matlines(ind3,res.block.dyn2$meanC[ind3,3],type="l",lty=3,lwd=1.5,col=par.col[2])
matlines(ind4,res.block.dyn4$meanC[ind4,3],type="l",lty=4,lwd=1.5,col=par.col[2])
matlines(ind2,res.block.dyn3$meanC[ind2,3],type="l",lty=5,lwd=1.5,col=par.col[2])
dev.off()

rm(ind,ind1,ind2,ind3,ind4)
rm(res.mcmc.dyn,res.block.dyn2,res.block.dyn3,res.block.dyn4,res.block.dyn1)


 #****************** Boxplot of the path functionals *****************************#

   
save.path.data = "folder.to.save.the.simulated.dataset"

N=1000; b.block=c(50,100,200)


save.dir.res.block =  paste(save.path.data,"Block/N=",N,"/",sep="")
save.dir.res.mcmc  =  paste(save.path.data,"MCMC/",sep="")
save.dir.res.smc   =  paste(save.path.data,"SMC/N=",N,"/",sep="")

 # time points
t1 = 10000; t2 = 20000; t3=30000

for(i in 1:50)
    {
    print(i)
    setwd(save.dir.res.block)
    res.block50 =   summary(readRDS(paste("data1_block_N=",N,"_b=",50,"_MC",i,".rds",sep="")))[4]  # meanC
    res.block100 =  summary(readRDS(paste("data1_block_N=",N,"_b=",100,"_MC",i,".rds",sep="")))[4]  # meanC
    res.block200 =  summary(readRDS(paste("data1_block_N=",N,"_b=",200,"_MC",i,".rds",sep="")))[4]  # meanC
    setwd(save.dir.res.smc)
    res.smc  = summary(readRDS(paste("data1_smc_N=",N,"_MC",i,".rds",sep="")))[4] # data1_smc_N=100_MCi

    #97.5%
    if(i==1){
rest1.Q2 = rest2.Q2 = rest3.Q2 = matrix(,1,4); rest1.mean = rest2.mean = rest3.mean = matrix(,1,4); 

rest1.Q2[1,1] = res.block50[[1]][t1,3];  rest2.Q2[1,1] = res.block50[[1]][t2,3]; rest3.Q2[1,1] = res.block50[[1]][t3,3];
rest1.Q2[1,2] = res.block100[[1]][t1,3]; rest2.Q2[1,2] = res.block100[[1]][t2,3];rest3.Q2[1,2] = res.block100[[1]][t3,3]
rest1.Q2[1,3] = res.block200[[1]][t1,3]; rest2.Q2[1,3] = res.block200[[1]][t2,3];rest3.Q2[1,3] = res.block200[[1]][t3,3]; 
rest1.Q2[1,4] = res.smc[[1]][t1,3];      rest2.Q2[1,4] = res.smc[[1]][t2,3]; rest3.Q2[1,4] = res.smc[[1]][t3,3];

rest1.mean[1,1] = res.block50[[1]][t1,1];  rest2.mean[1,1] = res.block50[[1]][t2,1]; rest3.mean[1,1] = res.block50[[1]][t3,1];
rest1.mean[1,2] = res.block100[[1]][t1,1]; rest2.mean[1,2] = res.block100[[1]][t2,1]; rest3.mean[1,2] = res.block100[[1]][t3,1]
rest1.mean[1,3] = res.block200[[1]][t1,1]; rest2.mean[1,3] = res.block200[[1]][t2,1]; rest3.mean[1,3] = res.block200[[1]][t3,1]; 
rest1.mean[1,4] = res.smc[[1]][t1,1];      rest2.mean[1,4] = res.smc[[1]][t2,1]; rest3.mean[1,4] = res.smc[[1]][t3,1];
            }
    if(i >1)
        {
rest1.Q2  = rbind(rest1.Q2,c(res.block50[[1]][t1,3],res.block100[[1]][t1,3],res.block200[[1]][t1,3],res.smc[[1]][t1,3]))
rest2.Q2  = rbind(rest2.Q2,c(res.block50[[1]][t2,3],res.block100[[1]][t2,3],res.block200[[1]][t2,3],res.smc[[1]][t2,3]))
rest3.Q2  = rbind(rest3.Q2,c(res.block50[[1]][t3,3],res.block100[[1]][t3,3],res.block200[[1]][t3,3],res.smc[[1]][t3,3]))

rest1.mean  = rbind(rest1.mean,c(res.block50[[1]][t1,1],res.block100[[1]][t1,1],res.block200[[1]][t1,1],res.smc[[1]][t1,1]))
rest2.mean  = rbind(rest2.mean,c(res.block50[[1]][t2,1],res.block100[[1]][t2,1],res.block200[[1]][t2,1],res.smc[[1]][t2,1]))
rest3.mean  = rbind(rest3.mean,c(res.block50[[1]][t3,1],res.block100[[1]][t3,1],res.block200[[1]][t3,1],res.smc[[1]][t3,1]))
    	}
rm(res.block50,res.block100,res.block200)
print(i)
}

setwd(save.dir.res.mcmc)
res.mcmc      = summary(readRDS("data1_MCMC.rds"))[4]
res.mcmc.Q2   = c(res.mcmc[[1]][t1,3],res.mcmc[[1]][t2,3],res.mcmc[[1]][t3,3])
res.mcmc.mean = c(res.mcmc[[1]][t1,1],res.mcmc[[1]][t2,1],res.mcmc[[1]][t3,1])


setwd(save.path.data)
pdf(paste("box1_Q2_N=",N,".pdf",sep=""))
main = paste("Time = ",t1,sep=""); par(cex.axis=1.5)
boxplot(rest1.Q2, ylim=NULL, names = c("block 50","block 100","block 200","PF"))
abline(h = res.mcmc.Q2[1], col = "red") 
dev.off()

pdf(paste("box2_Q2_N=",N,".pdf",sep=""))
main = paste("Time = ",t2,sep=""); par(cex.axis=1.5) 
boxplot(rest2.Q2, ylim=NULL, names = c("block 50","block 100","block 200","PF"))
abline(h = res.mcmc.Q2[2], col = "red") 
dev.off()

pdf(paste("box3_Q2_N=",N,".pdf",sep=""))
main = paste("Time = ",t3,sep=""); par(cex.axis=1.5) 
boxplot(rest3.Q2, ylim=NULL, names = c("block 50","block 100","block 200","PF"))
abline(h = res.mcmc.Q2[3], col = "red") 
dev.off()

setwd(save.path.data)
pdf(paste("box1_mean_N=",N,".pdf",sep=""))
main = paste("Time = ",t1,sep=""); par(cex.axis=1.5) 
boxplot(rest1.mean, ylim=NULL, names = c("block 50","block 100","block 200","PF"))
abline(h = res.mcmc.mean[1], col = "red") 
dev.off()

pdf(paste("box2_mean_N=",N,".pdf",sep=""))
main = paste("Time = ",t2,sep=""); par(cex.axis=1.5) 
boxplot(rest2.mean, ylim=NULL, names = c("block 50","block 100","block 200","PF"))
abline(h = res.mcmc.mean[2], col = "red") 
dev.off()

pdf(paste("box3_mean_N=",N,".pdf",sep=""))
main = paste("Time = ",t3,sep=""); par(cex.axis=1.5) 
boxplot(rest3.mean, ylim=NULL, names = c("block 50","block 100","block 200","PF"))
abline(h = res.mcmc.mean[3], col = "red") 
dev.off()



#************** Distance of the SMC path approx. and Kalman distribution************************#

memory.size(max = FALSE)
require(dllmsmc)

save.path.data = "folder.to.save.the.simulated.dataset"

N.sim=c(100,1000); b.block=c(50,100,200)

for(k in 1:2)
    {
N = N.sim[k]
    for(kk in 1:3)
        {
b = b.block[kk]

print("b = "); print(b); print("")     
save.dir.res.block =  paste(save.path.data,"Block/N=",N,"/",sep="")
save.dir.res.mcmc  =  paste(save.path.data,"MCMC/",sep="")
save.dir.res.smc   =  paste(save.path.data,"SMC/N=",N,"/",sep="")

   
  # Blocking   (mean of MC.runs,  ks.dllm)
     
     
setwd(save.dir.res.block)
name.file= paste("data1_block_N=",N,"_b=",b,"_MC",1,".rds",sep="") #data1_block_N=100_b=50_MC1
files <- list.files(path = save.dir.res.block)
n.files <-  length(files)-3 # should be equal to MC runs (# of .rds files)

print("number of files = "); print(n.files); print("")
n.files = 50

temp            = readRDS(name.file)
n.list.elements = length(names(temp))

res.dyn.sum   = summary(temp)[2:4] # "meanT"    "meanSQ"   "meanC"
res.dyn.sumSQ = lapply(res.dyn.sum,function(x) x^2)

setwd(save.dir.res.mcmc)
res.mcmc   = readRDS("data1_MCMC.rds")
dist.dyn   = ks.dllm(temp,res.mcmc)[1:3]   #"meanT"   "meanSQ"  "meanC"  
dist.dyn   = lapply(dist.dyn,function(x) x/50)
rm(temp)

setwd(save.dir.res.block)
for(i in 2:n.files) 
    {
  temp     = readRDS(paste("data1_block_N=",N,"_b=",b,"_MC",i,".rds",sep=""))
  temp.dyn = summary(temp)[2:4]; temp.dist = ks.dllm(temp,res.mcmc)[1:3]    
  rm(temp)
  for (j in 1:3)
        {
    #res.dyn[[j]] <- res.dyn[[j]] + temp.dyn[[j]]/50     
    dist.dyn[[j]]      <- dist.dyn[[j]]      + temp.dist[[j]]/50 
    res.dyn.sum[[j]]   <- res.dyn.sum[[j]]   + temp.dyn[[j]]
    res.dyn.sumSQ[[j]] <- res.dyn.sumSQ[[j]] + temp.dyn[[j]]^2    
        }   
    rm(temp.dyn,temp.dist); print("MC simulation = "); print(i); print("")
    print("b=");print(b);print(N)
    }
rm(res.mcmc)

names(dist.dyn)  =  c("meanT","meanSQ","meanC")
setwd(save.path.data);
saveRDS(dist.dyn,paste("dist_MCruns_block_N=",N,"_b=",b,".rds",sep="")); rm(dist.dyn)
print(""); print("done")
 }
	}





