set.seed(10001)
eta1.true.MAT <- matrix(rt(1000*11,df=5),nrow=1000,ncol=11)
eta0.true.MAT <- matrix(rt(1000*11,df=5),nrow=1000,ncol=11)
omega.true.MAT <- matrix(rt(1000*10,df=5),nrow=1000,ncol=10)
sigma2eps.true <- c(1,1)
sigma2sg.true <- c(1,1)

N <- 1000; Ntest=1000; MCsize <- 10000

mopt1.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt2.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt3.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt4.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt5.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
a1.opt.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt1.mean.emplik.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt2.mean.emplik.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt3.mean.emplik.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt4.mean.emplik.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt5.mean.emplik.MAT <- matrix(0,nrow=1000,ncol=Ntest)
a1.opt.emplik.MAT <- matrix(0,nrow=1000,ncol=Ntest)
QL.mopt1.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
QL.mopt2.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
QL.mopt3.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
QL.mopt4.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
QL.mopt5.mean.MAT <- matrix(0,nrow=1000,ncol=Ntest)
QL.a1.opt.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt1.true.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt2.true.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt3.true.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt4.true.MAT <- matrix(0,nrow=1000,ncol=Ntest)
mopt5.true.MAT <- matrix(0,nrow=1000,ncol=Ntest)
a1.opt.true.MAT <- matrix(0,nrow=1000,ncol=Ntest)

EV.Proposed <- rep(0,1000)
EV.Proposed.Emplik <- rep(0,1000)
EV.QL <- rep(0,1000)

Time.Proposed <- rep(0,1000)
Time.QL <- rep(0,1000)

xi.ID.vec <- sample(1:3,size=(Ntest*1000),replace=TRUE, prob = c(0.25,0.6,0.1))
xi.vec <- rep(0,(Ntest*1000))
xi.vec[xi.ID.vec==1] <- rnorm(sum(xi.ID.vec==1), mean=-4, sd=1)
xi.vec[xi.ID.vec==2] <- rnorm(sum(xi.ID.vec==2), mean=-2, sd=1)
xi.vec[xi.ID.vec==3] <- rnorm(sum(xi.ID.vec==3), mean=22, sd=1)
xiMAT.test <- matrix(xi.vec,nrow=1000,ncol=Ntest)


for(it in 1:200){
  
  set.seed(it)
  
  ###Generate random coefficient#####
  eta1.true <- eta1.true.MAT[it,]
  eta0.true <- eta0.true.MAT[it,]
  eta1MC.mock <- matrix(rep(eta1.true,2),nrow=2, ncol=11,byrow=T)
  eta0MC.mock <- matrix(rep(eta0.true,2),nrow=2, ncol=11,byrow=T)
  omegaMC.mock <- matrix(rep(omega.true.MAT[it,],2),nrow=2,ncol=10,byrow=T)
  
  ###Generate training data#####
  A1 <- sample(1:5,size = N, replace = TRUE)
  A2 <- rbinom(N,1,prob=0.5)
  A2id1 <- which(A2==1)
  A2id0 <- which(A2==0)
  
  X1 <- rt(N,df=10)
  STAGE1design <- cbind((A1==1), (A1==1)*X1, (A1==2), (A1==2)*X1, (A1==3), (A1==3)*X1, (A1==4), (A1==4)*X1, (A1==5), (A1==5)*X1)
  mix.err.ID <- sample(1:3,size=N,replace=TRUE, prob = c(0.25,0.6,0.1))
  X2.err <- rep(0,N)
  X2.err[mix.err.ID==1] <- rnorm(sum(mix.err.ID==1), mean=-4, sd=1)
  X2.err[mix.err.ID==2] <- rnorm(sum(mix.err.ID==2), mean=-2, sd=1)
  X2.err[mix.err.ID==3] <- rnorm(sum(mix.err.ID==3), mean=22, sd=1)
  X2 <- c(STAGE1design %*% omega.true.MAT[it,]) + X2.err
  STAGE2design <- cbind( apply(STAGE1design,2,function(r){ r*(A2==1)}),X2*(A2==1),apply(STAGE1design,2,function(r){ r*(A2==0)}),X2*(A2==0) )
  Y <- c(STAGE2design%*%c(eta1.true,eta0.true)) + rnorm(N)
  
  ####Generate test data
  set.seed(20000+it)
  X1.test <- rt(Ntest,df=10)
  
  
  for(i in 1:Ntest){
    
    zeta.ind.temp1 <- sample(1:3,size=MCsize,replace=TRUE, prob = c(0.25,0.6,0.1))
    zeta.temp1 <- rep(0,MCsize)
    zeta.temp1[zeta.ind.temp1==1] <- rnorm(sum(zeta.ind.temp1==1), mean=-4, sd=1)
    zeta.temp1[zeta.ind.temp1==2] <- rnorm(sum(zeta.ind.temp1==2), mean=-2, sd=1)
    zeta.temp1[zeta.ind.temp1==3] <- rnorm(sum(zeta.ind.temp1==3), mean=22, sd=1)
    EAbs1.temp <- mean(abs(eta1.true[1] - eta0.true[1] + (eta1.true[2]-eta0.true[2])*X1.test[i] + ( eta1.true[11]-eta0.true[11] )*( omega.true.MAT[it,1] + omega.true.MAT[it,2]*X1.test[i] + zeta.temp1 )))
    
    zeta.ind.temp2 <- sample(1:3,size=MCsize,replace=TRUE, prob = c(0.25,0.6,0.1))
    zeta.temp2 <- rep(0,MCsize)
    zeta.temp2[zeta.ind.temp2==1] <- rnorm(sum(zeta.ind.temp2==1), mean=-4, sd=1)
    zeta.temp2[zeta.ind.temp2==2] <- rnorm(sum(zeta.ind.temp2==2), mean=-2, sd=1)
    zeta.temp2[zeta.ind.temp2==3] <- rnorm(sum(zeta.ind.temp2==3), mean=22, sd=1)
    EAbs2.temp <- mean(abs(eta1.true[3] - eta0.true[3] + (eta1.true[4]-eta0.true[4])*X1.test[i] + ( eta1.true[11]-eta0.true[11] )*( omega.true.MAT[it,3] + omega.true.MAT[it,4]*X1.test[i] + zeta.temp2 )))
    
    zeta.ind.temp3 <- sample(1:3,size=MCsize,replace=TRUE, prob = c(0.25,0.6,0.1))
    zeta.temp3 <- rep(0,MCsize)
    zeta.temp3[zeta.ind.temp3==1] <- rnorm(sum(zeta.ind.temp3==1), mean=-4, sd=1)
    zeta.temp3[zeta.ind.temp3==2] <- rnorm(sum(zeta.ind.temp3==2), mean=-2, sd=1)
    zeta.temp3[zeta.ind.temp3==3] <- rnorm(sum(zeta.ind.temp3==3), mean=22, sd=1)
    EAbs3.temp <- mean(abs(eta1.true[5] - eta0.true[5] + (eta1.true[6]-eta0.true[6])*X1.test[i] + ( eta1.true[11]-eta0.true[11] )*( omega.true.MAT[it,5] + omega.true.MAT[it,6]*X1.test[i] + zeta.temp3 )))
    
    zeta.ind.temp4 <- sample(1:3,size=MCsize,replace=TRUE, prob = c(0.25,0.6,0.1))
    zeta.temp4 <- rep(0,MCsize)
    zeta.temp4[zeta.ind.temp4==1] <- rnorm(sum(zeta.ind.temp4==1), mean=-4, sd=1)
    zeta.temp4[zeta.ind.temp4==2] <- rnorm(sum(zeta.ind.temp4==2), mean=-2, sd=1)
    zeta.temp4[zeta.ind.temp4==3] <- rnorm(sum(zeta.ind.temp4==3), mean=22, sd=1)
    EAbs4.temp <- mean(abs(eta1.true[7] - eta0.true[7] + (eta1.true[8]-eta0.true[8])*X1.test[i] + ( eta1.true[11]-eta0.true[11] )*( omega.true.MAT[it,7] + omega.true.MAT[it,8]*X1.test[i] + zeta.temp4 )))
    
    zeta.ind.temp5 <- sample(1:3,size=MCsize,replace=TRUE, prob = c(0.25,0.6,0.1))
    zeta.temp5 <- rep(0,MCsize)
    zeta.temp5[zeta.ind.temp5==1] <- rnorm(sum(zeta.ind.temp5==1), mean=-4, sd=1)
    zeta.temp5[zeta.ind.temp5==2] <- rnorm(sum(zeta.ind.temp5==2), mean=-2, sd=1)
    zeta.temp5[zeta.ind.temp5==3] <- rnorm(sum(zeta.ind.temp5==3), mean=22, sd=1)
    EAbs5.temp <- mean(abs(eta1.true[9] - eta0.true[9] + (eta1.true[10]-eta0.true[10])*X1.test[i] + ( eta1.true[11]-eta0.true[11] )*( omega.true.MAT[it,9] + omega.true.MAT[it,10]*X1.test[i] + zeta.temp5 )))
    
    mopt1.true.MAT[it,i] <- 0.5*(eta1.true[1] + eta0.true[1] + (eta1.true[2]+eta0.true[2])*X1.test[i] + (eta1.true[11]+eta0.true[11])*( omega.true.MAT[it,1] + omega.true.MAT[it,2]*X1.test[i]) + EAbs1.temp )
    mopt2.true.MAT[it,i] <- 0.5*(eta1.true[3] + eta0.true[3] + (eta1.true[4]+eta0.true[4])*X1.test[i] + (eta1.true[11]+eta0.true[11])*( omega.true.MAT[it,3] + omega.true.MAT[it,4]*X1.test[i]) + EAbs2.temp )
    mopt3.true.MAT[it,i] <- 0.5*(eta1.true[5] + eta0.true[5] + (eta1.true[6]+eta0.true[6])*X1.test[i] + (eta1.true[11]+eta0.true[11])*( omega.true.MAT[it,5] + omega.true.MAT[it,6]*X1.test[i]) + EAbs3.temp )
    mopt4.true.MAT[it,i] <- 0.5*(eta1.true[7] + eta0.true[7] + (eta1.true[8]+eta0.true[8])*X1.test[i] + (eta1.true[11]+eta0.true[11])*( omega.true.MAT[it,7] + omega.true.MAT[it,8]*X1.test[i]) + EAbs4.temp )
    mopt5.true.MAT[it,i] <- 0.5*(eta1.true[9] + eta0.true[9] + (eta1.true[10]+eta0.true[10])*X1.test[i] + (eta1.true[11]+eta0.true[11])*( omega.true.MAT[it,9] + omega.true.MAT[it,10]*X1.test[i]) + EAbs5.temp )
    a1.opt.true.MAT[it,i] <- which.max(c(mopt1.true.MAT[it,i],mopt2.true.MAT[it,i],mopt3.true.MAT[it,i],mopt4.true.MAT[it,i],mopt5.true.MAT[it,i]))
    
  }
  
  
  
  ###Inference for proposed method#############################
  ###Distributed fork for faster computing#####################
  #############################################################
  start.Proposed <- Sys.time()
  LMobj <-lm(Y~STAGE2design-1)
  eta.hat <- coef(LMobj)
  RSS <- sum(summary(LMobj)$res^2)
  X2LMobj <- lm(X2~STAGE1design-1)
  omega.hat <- coef(X2LMobj)
  RSS.X2 <- sum(summary(X2LMobj)$res^2)
  
  Eta.MC <- DistributedComp.ctsX.sim2b(Y = Y, x2 = X2, stage1design = STAGE1design, stage2design = STAGE2design, etaHat = eta.hat, omegaHat = omega.hat, RSS.Y = RSS, RSS.X2 = RSS.X2, YorX2 = 1)
  Eta.samples <- Eta.MC$Tsamples
  Omega.MC <- DistributedComp.ctsX.sim2b(Y = Y, x2 = X2, stage1design = STAGE1design, stage2design = STAGE2design, etaHat = eta.hat, omegaHat = omega.hat, RSS.Y = RSS, RSS.X2 = RSS.X2, YorX2 = 0)
  sigma2g.samples <- Omega.MC$sigma2g.samples
  Omega.samples <- Omega.MC$Tsamples.X2
  
  
  ####Proposed method: estimate optimal treatment###########
  for(i in 1:Ntest){
    
    OptObj <- OptTrt.sim2b(etaHMC = Eta.samples, omegaHMC = Omega.samples, sigmasqHMC = sigma2g.samples, x1 = X1.test[i])
    mopt1.mean.MAT[it,i] <- OptObj$m.opt.stage1[1]
    mopt2.mean.MAT[it,i] <- OptObj$m.opt.stage1[2]
    mopt3.mean.MAT[it,i] <- OptObj$m.opt.stage1[3]
    mopt4.mean.MAT[it,i] <- OptObj$m.opt.stage1[4]
    mopt5.mean.MAT[it,i] <- OptObj$m.opt.stage1[5]
    a1.opt.MAT[it,i] <- getmode(OptObj$optA1)
    
  }
  end.Proposed <- Sys.time()
  Time.Proposed[it] <- as.numeric(difftime(end.Proposed,start.Proposed,units = "secs"))
  
  STAGE1design.proposed <- cbind((a1.opt.MAT[it,]==1), (a1.opt.MAT[it,]==1)*X1.test, (a1.opt.MAT[it,]==2), (a1.opt.MAT[it,]==2)*X1.test, (a1.opt.MAT[it,]==3), (a1.opt.MAT[it,]==3)*X1.test, (a1.opt.MAT[it,]==4), (a1.opt.MAT[it,]==4)*X1.test, (a1.opt.MAT[it,]==5), (a1.opt.MAT[it,]==5)*X1.test)
  X2.proposed <- c(STAGE1design.proposed %*% omega.true.MAT[it,]) + xiMAT.test[it,]
  a2.test.proposed <- as.numeric(c(cbind(STAGE1design.proposed,X2.proposed)%*%eta1.true) >= c(cbind(STAGE1design.proposed,X2.proposed)%*%eta0.true))
  STAGE2design.proposed <- cbind( apply(STAGE1design.proposed,2,function(r){ r*(a2.test.proposed==1)}),X2.proposed*(a2.test.proposed==1),apply(STAGE1design.proposed,2,function(r){ r*(a2.test.proposed==0)}),X2.proposed*(a2.test.proposed==0) )
  EV.Proposed[it] <- mean(c(STAGE2design.proposed%*%c(eta1.true,eta0.true)))
  
  Eta.emplik.MC <- DistributedComp.ctsX.sim2b.Emplik.VB( y = Y, x2 = X2, stage1design = STAGE1design, stage2design = STAGE2design, etaHat = eta.hat, omegaHat = omega.hat, YorX2 = 1, nsamp = 10000)
  Omega.emplik.MC <- DistributedComp.ctsX.sim2b.Emplik.VB( y = Y, x2 = X2, stage1design = STAGE1design, stage2design = STAGE2design, etaHat = eta.hat, omegaHat = omega.hat, YorX2 = 2, nsamp = 10000)
  #WeightMatrix <- matrix(0,nrow=nrow(Omega.emplik.MC$Tsamples),ncol=N)
  ZetaMatrix <- matrix(0,nrow=nrow(Omega.emplik.MC$Tsamples),ncol=N)
  for(j in 1:nrow(Omega.emplik.MC$Tsamples)){
    
    #WeightMatrix[j,] <- evalGel(g=FSLR, x = cbind(STAGE1design,X2), tet0 = Omega.emplik.MC$Tsamples[j,], type="EL", optlam = "optim")$pt
    ZetaMatrix[j,] <- X2 - c(STAGE1design%*%Omega.emplik.MC$Tsamples[j,])
    
  }
  for(i in 1:Ntest){
    
    OptObj.emplik <- OptTrt.sim2b.emplik(omegaHMC = Omega.emplik.MC$Tsamples[1:10000,], etaHMC = Eta.emplik.MC$Tsamples[1:10000,], x1 = X1.test[i], XiMatrix = ZetaMatrix[1:10000,])
    mopt1.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[1]
    mopt2.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[2]
    mopt3.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[3]
    mopt4.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[4]
    mopt5.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[5]
    a1.opt.emplik.MAT[it,i] <- getmode(OptObj.emplik$optA1)
    
  }
  
  STAGE1design.proposed.emplik <- cbind((a1.opt.emplik.MAT[it,]==1), (a1.opt.emplik.MAT[it,]==1)*X1.test, (a1.opt.emplik.MAT[it,]==2), (a1.opt.emplik.MAT[it,]==2)*X1.test, (a1.opt.emplik.MAT[it,]==3), (a1.opt.emplik.MAT[it,]==3)*X1.test, (a1.opt.emplik.MAT[it,]==4), (a1.opt.emplik.MAT[it,]==4)*X1.test, (a1.opt.emplik.MAT[it,]==5), (a1.opt.emplik.MAT[it,]==5)*X1.test)
  X2.proposed.emplik <- c(STAGE1design.proposed.emplik %*% omega.true.MAT[it,]) + xiMAT.test[it,]
  a2.test.proposed.emplik <- as.numeric(c(cbind(STAGE1design.proposed.emplik,X2.proposed.emplik)%*%eta1.true) >= c(cbind(STAGE1design.proposed.emplik,X2.proposed.emplik)%*%eta0.true))
  STAGE2design.proposed.emplik <- cbind( apply(STAGE1design.proposed.emplik,2,function(r){ r*(a2.test.proposed.emplik==1)}),X2.proposed.emplik*(a2.test.proposed.emplik==1),apply(STAGE1design.proposed.emplik,2,function(r){ r*(a2.test.proposed.emplik==0)}),X2.proposed.emplik*(a2.test.proposed.emplik==0) )
  EV.Proposed.Emplik[it] <- mean(c(STAGE2design.proposed.emplik%*%c(eta1.true,eta0.true)))
  
  
  ###Q-Learning#####
  dataset.QL <- cbind(unname(X1),unname(X2),A1,A2, Y)
  colnames(dataset.QL) <- c("o1","o2","a1","a2","y")
  start.QL <- Sys.time()
  QLobj <- ql.glm.ctsX.sim2b(dataset=list(data=dataset.QL), o1.test = X1.test, stage2design=STAGE2design, stage1design=STAGE1design)
  QL.mopt1.mean.MAT[it,] <- QLobj$mu1.hat[,1]
  QL.mopt2.mean.MAT[it,] <- QLobj$mu1.hat[,2]
  QL.mopt3.mean.MAT[it,] <- QLobj$mu1.hat[,3]
  QL.mopt4.mean.MAT[it,] <- QLobj$mu1.hat[,4]
  QL.mopt5.mean.MAT[it,] <- QLobj$mu1.hat[,5]
  QL.a1.opt.MAT[it,] <- QLobj$a1opt
  end.QL <- Sys.time()
  Time.QL[it] <- as.numeric(difftime(end.QL,start.QL,units = "secs"))
  STAGE1design.QL <- cbind((QL.a1.opt.MAT[it,]==1), (QL.a1.opt.MAT[it,]==1)*X1.test, (QL.a1.opt.MAT[it,]==2), (QL.a1.opt.MAT[it,]==2)*X1.test, (QL.a1.opt.MAT[it,]==3), (QL.a1.opt.MAT[it,]==3)*X1.test, (QL.a1.opt.MAT[it,]==4), (QL.a1.opt.MAT[it,]==4)*X1.test, (QL.a1.opt.MAT[it,]==5), (QL.a1.opt.MAT[it,]==5)*X1.test)
  X2.QL <- c(STAGE1design.QL %*% omega.true.MAT[it,]) + xiMAT.test[it,]
  a2.QL <- as.numeric(c(cbind(STAGE1design.QL,X2.QL)%*%eta1.true) >= c(cbind(STAGE1design.QL,X2.QL)%*%eta0.true))
  STAGE2design.QL <- cbind( apply(STAGE1design.QL,2,function(r){ r*(a2.QL==1)}),X2.QL*(a2.QL==1),apply(STAGE1design.QL,2,function(r){ r*(a2.QL==0)}),X2.QL*(a2.QL==0) )
  EV.QL[it] <- mean(c(STAGE2design.QL%*%c(eta1.true,eta0.true)))
  
  
  rm(OptObj)
  rm(OptObj.emplik)
  #rm(OptObj.true)
  rm(QLobj)
  #rm(BMLobj)
  gc(); gc()
  if(it%%10==0){
    
    cat("Completed it=", it, "\n")
    save.image("Simulation1EmpLikDynTrtRevisionN1000v1.RData")
    
  }
  
}
