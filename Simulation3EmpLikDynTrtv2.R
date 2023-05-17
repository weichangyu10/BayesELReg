#Simulate true coefficient values for all 1000 repetitions
ETA.func <- function(x){
  
  1/(2^3 *pi)*x^(2)*log(x)
  
}

m1star <- function(x1,a1,x2){
  
  cos(pi*x2*0.5 + 0.7*pi*x1*a1) + 0.1*a1
  
}

m0star <- function(x1,a1,x2){
  
  -cos(pi*x2*0.75 - pi*x1) + 0.25*a1
  
}

gstar <- function(x1, a1){
  
  a1*(2*expit(6*sin(x1^2)) - 1)
  
}

N <- 2000; N.test <- 100; N.EV <- 200

#Initialise storage objects
mopt1.mean.MAT <- matrix(0,nrow=1000,ncol=N.test)
moptNeg1.mean.MAT <- matrix(0,nrow=1000,ncol=N.test)
mopt1.mean.emplik.MAT <- matrix(0,nrow=1000,ncol=N.test)
moptNeg1.mean.emplik.MAT <- matrix(0,nrow=1000,ncol=N.test)
mopt1.true.MAT <- matrix(0,nrow=1000,ncol=N.test)
moptNeg1.true.MAT <- matrix(0,nrow=1000,ncol=N.test)
a1.opt.MAT <- matrix(0,nrow=1000,ncol=N.test)
a1.opt.emplik.MAT <- matrix(0,nrow=1000,ncol=N.test)
a1.opt.true.MAT <- matrix(0,nrow=1000,ncol=N.test)
QL.a1.opt.MAT <- matrix(0,nrow=1000,ncol=N.test)
QL.mu1hat.opt.MAT <- matrix(0,nrow=1000,ncol=(2*N.test))
BML.a1.opt.MAT <- matrix(0,nrow=1000,ncol=N.test)
BML.mu1hat.opt.MAT <- matrix(0,nrow=1000,ncol=(2*N.test))

Time.Proposed <- rep(0,1000)
Time.QL <- rep(0,1000)
Time.BML <- rep(0,1000)

cores=detectCores()

EV.Proposed <- rep(0,1000)
EV.Proposed.Emplik <- rep(0,1000)
EV.QL <- rep(0,1000)
EV.BML <- rep(0,1000)

xi.ID.vec <- sample(1:3,size=(N.test*1000),replace=TRUE, prob = c(0.25,0.6,0.1))
xi.vec <- rep(0,(N.test*1000))
xi.vec[xi.ID.vec==1] <- rnorm(sum(xi.ID.vec==1), mean=-0.4, sd=0.1)
xi.vec[xi.ID.vec==2] <- rnorm(sum(xi.ID.vec==2), mean=-0.2, sd=0.1)
xi.vec[xi.ID.vec==3] <- rnorm(sum(xi.ID.vec==3), mean=2.2, sd=0.1)
xi.test.MAT <- matrix(xi.vec,nrow=1000,ncol=N.test)
epsilon.test.MAT <- matrix(rnorm((N.test*1000),mean=0,sd=0.25),nrow=1000,ncol=N.test)

for(it in 1:200){
  
  set.seed(it)
  knots.choice <- 5
  knots.choice.proposed <- 5
  ##Generate training data
  x1 <- runif(N, min = -1, max= 1)
  a1 <- 2*rbinom(N,size=1,prob=0.5) - 1
  x2.err.gp <- sample(1:3,size=N,replace=TRUE, prob = c(0.25,0.6,0.1))
  x2.err <- rep(0,N)
  x2.err[x2.err.gp==1] <- rnorm(sum(x2.err.gp==1), mean=-0.4, sd=0.1)
  x2.err[x2.err.gp==2] <- rnorm(sum(x2.err.gp==2), mean=-0.2, sd=0.1)
  x2.err[x2.err.gp==3] <- rnorm(sum(x2.err.gp==3), mean=2.2, sd=0.1)
  x2 <- gstar(x1=x1,a1=a1) + x2.err
  a2 <- sort(rbinom(N,size=1,prob=0.5),decreasing = TRUE)
  n1 <- sum(a2); n0 <- N - n1
  y <- rep(0,N)
  y[1:n1] <- m1star(x1=x1[1:n1],a1=a1[1:n1],x2=x2[1:n1]) + rnorm(n1,sd=0.25)
  y[-(1:n1)] <- m0star(x1=x1[-(1:n1)],a1=a1[-(1:n1)],x2=x2[-(1:n1)]) + rnorm(n0,sd=0.25)
  stage1knots <- matrix(quantile(x1,seq(0,0.975,length.out=knots.choice)), ncol=1)
  stage1knots.proposed <- matrix(quantile(x1,seq(0,0.975,length.out=knots.choice.proposed)), ncol=1)
  stage2knots <- unname(as.matrix(expand.grid(stage1knots,quantile(x2,seq(0,0.975,length.out=knots.choice)))))
  stage2knots.proposed <- unname(as.matrix(expand.grid(stage1knots.proposed,quantile(x2,seq(0,0.975,length.out=knots.choice.proposed)))))
  BASIS1obj <- ns(x1,knots=stage1knots)
  BASIS1obj.proposed <- ns(x1,knots=stage1knots.proposed)
  
  stage2.basis <- basis.tps(x=cbind(x1,x2),knots=stage2knots, intercept=TRUE)
  stage2.basis.proposed <- basis.tps(x=cbind(x1,x2),knots=stage2knots.proposed, intercept=TRUE)
  stage2.predictors <- cbind(apply(stage2.basis,2,function(s){s*(a1==-1)*(a2==1)}),apply(stage2.basis,2,function(s){s*(a1==1)*(a2==1)}),apply(stage2.basis,2,function(s){s*(a1==-1)*(a2==0)}),apply(stage2.basis,2,function(s){s*(a1==1)*(a2==0)}))
  stage2.predictors.proposed <- cbind(apply(stage2.basis.proposed,2,function(s){s*(a1==-1)*(a2==1)}),apply(stage2.basis.proposed,2,function(s){s*(a1==1)*(a2==1)}),apply(stage2.basis.proposed,2,function(s){s*(a1==-1)*(a2==0)}),apply(stage2.basis.proposed,2,function(s){s*(a1==1)*(a2==0)}))
  stage1.basis <- matrix(c(ns(x1,knots=stage1knots)),nrow=N,ncol=(knots.choice+1) )
  stage1.basis.proposed <- matrix(c(ns(x1,knots=stage1knots.proposed)),nrow=N,ncol=(knots.choice.proposed+1) )
  stage1.predictors <- cbind(apply(stage1.basis,2,function(s){s*(a1==-1)}),apply(stage1.basis,2,function(s){s*(a1==1)}))
  stage1.predictors.proposed <- cbind(apply(stage1.basis.proposed,2,function(s){s*(a1==-1)}),apply(stage1.basis.proposed,2,function(s){s*(a1==1)}))
  
  
  ####Generate test data######
  x1.test <- runif(N.test, min = -1, max= 1)
  stage1.test.basis <- matrix(c(predict(BASIS1obj, x1.test)),nrow=N.test,ncol=(knots.choice+1))
  stage1.test.basis.proposed <- matrix(c(predict(BASIS1obj.proposed, x1.test)),nrow=N.test,ncol=(knots.choice.proposed+1))
  
  #Fit proposed method
  start.Proposed <- Sys.time()
  eta.hat <- psolve(t(stage2.predictors.proposed)%*%stage2.predictors.proposed,t(stage2.predictors.proposed)%*%y)
  RSS.Y <- sum((y-c(stage2.predictors.proposed%*%eta.hat))^2)
  omega.hat <- solve(t(stage1.predictors.proposed)%*%stage1.predictors.proposed)%*%t(stage1.predictors.proposed)%*%x2
  RSS.X2 <- sum((x2-c(stage1.predictors.proposed%*%omega.hat))^2)
  
  Eta.samples <- DistributedComp.ctsX.sim4(Y = y, x2 =  x2, stage1design = stage1.predictors.proposed, stage2design = stage2.predictors.proposed, etaHat = eta.hat, omegaHat = omega.hat, RSS.Y = RSS.Y, RSS.X2 = RSS.X2, YorX2 = 1)
  Omega.samples <- DistributedComp.ctsX.sim4(Y = y, x2 =  x2, stage1design = stage1.predictors.proposed, stage2design = stage2.predictors.proposed, etaHat = eta.hat, omegaHat = omega.hat, RSS.Y = RSS.Y, RSS.X2 = RSS.X2, YorX2 = 2)
  
  finalMatrix <- matrix(0,nrow=3,ncol=N.test)
  for(i in 1:N.test){
    
    OptObj <- OptTrt.ctsX.sim4(eta1MC = Eta.samples$Tsamples[,1:(length(eta.hat)/2)], eta0MC = Eta.samples$Tsamples[,(length(eta.hat)/2+1):length(eta.hat)], omegaMC = Omega.samples$Tsamples.X2, sigma2gMC = Omega.samples$sigma2g.samples, KnotsMatrix = stage2knots.proposed, x1.test.point = x1.test[i], x1.basis = stage1.test.basis.proposed[i,] )
    finalMatrix[,i] <- c(OptObj$m1opt.s1,OptObj$m1opt.sNeg1,OptObj$optA1)
    
  }
  #setup parallel backend to use many processors
  # cl <- makeCluster(cores[1]-1) #not to overload your computer
  # registerDoParallel(cl)
  # finalMatrix <- foreach(i=1:N.test, .combine=cbind) %dopar% {
  #   
  #   library(npreg)
  #   OptObj <- OptTrt.ctsX.sim4(eta1MC = Eta.samples$Tsamples[,1:(length(eta.hat)/2)], eta0MC = Eta.samples$Tsamples[,(length(eta.hat)/2+1):length(eta.hat)], omegaMC = Omega.samples$Tsamples.X2, sigma2gMC = Omega.samples$sigma2g.samples, KnotsMatrix = stage2knots.proposed, x1.test.point = x1.test[i], x1.basis = stage1.test.basis.proposed[i,] )
  #   c(OptObj$m1opt.s1,OptObj$m1opt.sNeg1,OptObj$optA1)
  # 
  # }
  # stopCluster(cl)
  end.Proposed <- Sys.time()
  Time.Proposed[it] <- as.numeric(difftime(end.Proposed,start.Proposed,units = "secs"))
  mopt1.mean.MAT[it,] <-  finalMatrix[1,]
  moptNeg1.mean.MAT[it,] <- finalMatrix[2,]
  a1.opt.MAT[it,] <- finalMatrix[3,]
  q2 <- ncol(Eta.samples$Tsamples)
  
  xi.test <- xi.test.MAT[it,]
  epsilon.test <- epsilon.test.MAT[it,]
  
  x2.test.proposed <- gstar(x1=x1.test,a1=a1.opt.MAT[it,]) + xi.test
  stage2.test.basis.proposed <- basis.tps( x=cbind(x1.test,x2.test.proposed),  knots = stage2knots.proposed, intercept = TRUE  )
  stage2.test.basis.a1.proposed <- cbind(apply(stage2.test.basis.proposed,2,function(s){ s*(a1.opt.MAT[it,]==-1)  }), apply(stage2.test.basis.proposed,2,function(s){ s*(a1.opt.MAT[it,]==1)  }) )
  a2.test.proposed <- as.numeric(apply((stage2.test.basis.a1.proposed %*% t(Eta.samples$Tsamples[,1:(q2/2)])) >= (stage2.test.basis.a1.proposed %*% t(Eta.samples$Tsamples[,(q2/2+1):q2])),1,getmode))
  y.test.proposed <- rep(NA,N.test)
  y.test.proposed[a2.test.proposed==1] <- m1star(x1=x1.test[a2.test.proposed==1],a1=a1.opt.MAT[it,a2.test.proposed==1],x2=x2.test.proposed[a2.test.proposed==1]) + epsilon.test[a2.test.proposed==1]
  y.test.proposed[a2.test.proposed==0] <- m0star(x1=x1.test[a2.test.proposed==0],a1=a1.opt.MAT[it,a2.test.proposed==0],x2=x2.test.proposed[a2.test.proposed==0]) + epsilon.test[a2.test.proposed==0]
  EV.Proposed[it] <- mean(y.test.proposed)
    # cl <- makeCluster(cores[1]-1) #not to overload your computer
    # registerDoParallel(cl)
    # finalMatrix.train <- foreach(i=1:N, .combine=cbind) %dopar% {
    #   
    #   library(npreg)
    #   OptObj <- OptTrt.ctsX.sim4(eta1MC = Eta.samples$Tsamples[,1:(length(eta.hat)/2)], eta0MC = Eta.samples$Tsamples[,(length(eta.hat)/2+1):length(eta.hat)], omegaMC = Omega.samples$Tsamples.X2, sigma2gMC = Omega.samples$sigma2g.samples, KnotsMatrix = stage2knots.proposed, x1.test.point = x1[i], x1.basis = stage1.basis.proposed[i,] )
    #   c(OptObj$m1opt.s1,OptObj$m1opt.sNeg1,OptObj$optA1)
    #   
    # }
    # stopCluster(cl)
  

  #Empirical likelihood
  Omega.Emplik.samples <- DistributedComp.ctsX.sim3.Emplik.VB(y = y, x2 = unname(x2), stage1design = stage1.predictors.proposed, stage2design = stage2.predictors.proposed, etaHat = eta.hat, omegaHat = omega.hat, YorX2 = 2, nsamp = 5000, VBiter = 5000)
  ZetaMatrix <- matrix(0,nrow=nrow(Omega.Emplik.samples$Tsamples),ncol=N)
  for(j in 1:nrow(Omega.Emplik.samples$Tsamples)){
    
    ZetaMatrix[j,] <- x2 - c(stage1.predictors%*%Omega.Emplik.samples$Tsamples[j,])
    
  }
  
  
  for(i in 1:N.test){
    
    OptObj.emplik <- OptTrt.sim4.emplik(eta1MC = Eta.samples$Tsamples[,1:(ncol(Eta.samples$Tsamples)/2)], eta0MC = Eta.samples$Tsamples[,(ncol(Eta.samples$Tsamples)/2+1):ncol(Eta.samples$Tsamples)], omegaMC = Omega.Emplik.samples$Tsamples, KnotsMatrix = stage2knots, x1.test.point = x1.test[i], x1.basis = stage1.test.basis[i,], XiMatrix = ZetaMatrix)
    mopt1.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[1]
    moptNeg1.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[2]
    #mopt3.mean.emplik.MAT[it,i] <- OptObj.emplik$m.opt.stage1[3]
    a1.opt.emplik.MAT[it,i] <- getmode(OptObj.emplik$optA1)
    rm(OptObj.emplik)
    gc();gc();gc()
  }
  a1.opt.emplik.MAT[it,(a1.opt.emplik.MAT[it,]==1)] <- -1
  a1.opt.emplik.MAT[it,(a1.opt.emplik.MAT[it,]==2)] <- 1
  
  x2.test.proposed.emplik <- gstar(x1=x1.test,a1=a1.opt.emplik.MAT[it,]) + xi.test
  stage2.test.basis.proposed.emplik <- basis.tps( x=cbind(x1.test,x2.test.proposed.emplik),  knots = stage2knots.proposed, intercept = TRUE  )
  stage2.test.basis.a1.proposed.emplik <- cbind(apply(stage2.test.basis.proposed.emplik,2,function(s){ s*(a1.opt.emplik.MAT[it,]==-1)  }), apply(stage2.test.basis.proposed.emplik,2,function(s){ s*(a1.opt.emplik.MAT[it,]==1)  }) )
  a2.test.proposed.emplik <- as.numeric(apply((stage2.test.basis.a1.proposed.emplik %*% t(Eta.samples$Tsamples[,1:(q2/2)])) >= (stage2.test.basis.a1.proposed.emplik %*% t(Eta.samples$Tsamples[,(q2/2+1):q2])),1,getmode))
  y.test.proposed.emplik <- rep(NA,N.test)
  y.test.proposed.emplik[a2.test.proposed.emplik==1] <- m1star(x1=x1.test[a2.test.proposed.emplik==1],a1=a1.opt.emplik.MAT[it,a2.test.proposed.emplik==1],x2=x2.test.proposed.emplik[a2.test.proposed.emplik==1]) + epsilon.test[a2.test.proposed.emplik==1]
  y.test.proposed.emplik[a2.test.proposed.emplik==0] <- m0star(x1=x1.test[a2.test.proposed.emplik==0],a1=a1.opt.emplik.MAT[it,a2.test.proposed.emplik==0],x2=x2.test.proposed.emplik[a2.test.proposed.emplik==0]) + epsilon.test[a2.test.proposed.emplik==0]
  EV.Proposed.Emplik[it] <- mean(y.test.proposed.emplik)
  
  #True value
  x2.draw.true.s1 <- matrix(rep(gstar(x1=x1.test,a1=rep(1,N.test)),1000) + rnorm(N.test*1000,sd=0.25),nrow=N.test,ncol=1000)
  x1mass.true <- rep(x1.test,1000)
  m1starmass.s1 <- m1star(x1mass.true, rep(1,length(x1mass.true)),  c(x2.draw.true.s1))
  m0starmass.s1 <- m0star(x1mass.true, rep(1,length(x1mass.true)),  c(x2.draw.true.s1))
  max.star.mass.s1 <- matrix(0.5*( m1starmass.s1 + m0starmass.s1 + abs(m1starmass.s1 - m0starmass.s1) ), nrow=N.test,ncol=1000)
  EYopt.s1.true <- rowMeans(max.star.mass.s1)
  mopt1.true.MAT[it,] <- EYopt.s1.true
  
  x2.draw.true.sNeg1 <- matrix(rep(gstar(x1=x1.test,a1=rep(-1,N.test)),1000) + rnorm(N.test*1000,sd=0.25),nrow=N.test,ncol=1000)
  m1starmass.sNeg1 <- m1star(x1mass.true, rep(-1,length(x1mass.true)),  c(x2.draw.true.sNeg1))
  m0starmass.sNeg1 <- m0star(x1mass.true, rep(-1,length(x1mass.true)),  c(x2.draw.true.sNeg1))
  max.star.mass.sNeg1 <- matrix(pmax(m1starmass.sNeg1,m0starmass.sNeg1), nrow=N.test,ncol=1000)
  EYopt.sNeg1.true <- rowMeans(max.star.mass.sNeg1)
  moptNeg1.true.MAT[it,] <- EYopt.sNeg1.true
  a1.opt.true.MAT[it,] <- sign(EYopt.s1.true - EYopt.sNeg1.true)
  
  
  #Q-learning
  dataset.QL <- cbind(unname(x1),unname(x2),a1,a2, y)
  colnames(dataset.QL) <- c("o1","o2","a1","a2","y")
  start.QL <- Sys.time()
  QLobj <- ql.sim4(dataset = list(data=dataset.QL, stage1basis=stage1.basis, stage2basis = stage2.basis), o1.test = x1.test, stage1.test.basis = stage1.test.basis)
  q1 <- length(QLobj$beta1.hat)
  QL.a1.opt.MAT[it,] <- QLobj$a1.opt
  QL.mu1hat.opt.MAT[it,] <- QLobj$mu1.hat
  end.QL <- Sys.time()
  Time.QL[it] <- as.numeric(difftime(end.QL,start.QL,units = "secs"))
  
  x2.test.QL <- gstar(x1=x1.test,a1=QLobj$a1.opt) + xi.test
  stage2.test.basis.QL <- basis.tps( x=cbind(x1.test,x2.test.QL),  knots = stage2knots, intercept = TRUE  )
  stage2.test.basis.a1.QL <- cbind(apply(stage2.test.basis.QL,2,function(s){ s*(QLobj$a1.opt==-1)  }), apply(stage2.test.basis.QL,2,function(s){ s*(QLobj$a1.opt==1)  }) )
  a2.test.QL <- as.numeric( (stage2.test.basis.a1.QL %*% QLobj$beta2.hat[1:(q2/2)]) >=  (stage2.test.basis.a1.QL %*% QLobj$beta2.hat[(q2/2+1):q2])  )
  y.test.QL <- rep(NA,N.test)
  y.test.QL[a2.test.QL==1] <- m1star(x1=x1.test[a2.test.QL==1],a1=QL.a1.opt.MAT[it,a2.test.QL==1],x2=x2.test.QL[a2.test.QL==1]) + epsilon.test[a2.test.QL==1]
  y.test.QL[a2.test.QL==0] <- m0star(x1=x1.test[a2.test.QL==0],a1=QL.a1.opt.MAT[it,a2.test.QL==0],x2=x2.test.QL[a2.test.QL==0]) + epsilon.test[a2.test.QL==0]
  EV.QL[it] <- mean(y.test.QL)
  
  #BML
  
  
  
  rm(x2.draw.true.s1)
  rm(x1mass.true)
  rm(m1starmass.s1)
  rm(m0starmass.s1)
  rm(m1starmass.sNeg1)
  rm(m0starmass.sNeg1)
  rm(max.star.mass.s1)
  rm(max.star.mass.sNeg1)
  rm(stage2.basis)
  rm(stage2.predictors)
  rm(stage1.basis)
  rm(stage1.predictors)
  rm(Eta.samples)
  rm(Omega.samples)
  rm(QLobj)
  rm(finalMatrix)
  #rm(OptObj)
  gc();gc();
  
  if(it%%2==0){
    
    save.image("Simulation3EmpLikDynTrtRevisionN2000Newv1.RData")
    
  }
  
}

