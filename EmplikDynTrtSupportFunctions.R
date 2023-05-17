#library(mvtnorm)
#library(LaplacesDemon)
library(elhmc)
library(gmm)
library(linpk)

OptTrt.sim2b <- function(etaHMC, omegaHMC, sigmasqHMC, x1){
  
  B <- nrow(etaHMC)
  eta1HMC <- etaHMC[,1:(ncol(etaHMC)/2)]
  eta0HMC <- etaHMC[,(ncol(etaHMC)/2+1):(ncol(etaHMC))]
  A1.size <- ((ncol(etaHMC)-2)/4)
  eta13 <- eta1HMC[,ncol(eta1HMC)]
  eta03 <- eta0HMC[,ncol(eta0HMC)]
  EYoptMatrix <- matrix(0,nrow=B,ncol=A1.size)
  for(s in 1:A1.size){
    
    eta11.temp <- eta1HMC[,((s-1)*2+1)]
    eta01.temp <- eta0HMC[,((s-1)*2+1)]
    eta12.temp <- eta1HMC[,(s*2)]
    eta02.temp <- eta0HMC[,(s*2)]
    omega1.temp <- omegaHMC[,((s-1)*2+1)]
    omega2.temp <- omegaHMC[,(s*2)]
    m.u <- (eta11.temp - eta01.temp) + (eta13 - eta03)*omega1.temp + ((eta12.temp - eta02.temp) + (eta13 - eta03)*omega2.temp)*x1
    EYoptMatrix[,s] <- 0.5*( (eta11.temp + eta01.temp) + (eta13 + eta03)*omega1.temp + ( (eta12.temp + eta02.temp) + omega2.temp*(eta13 + eta03))*x1 + sqrt(sigmasqHMC)*sqrt(2/pi) * exp( -m.u^2/(2*sigmasqHMC) ) + m.u*( 1- 2*pnorm(-m.u/sqrt(sigmasqHMC)) ) )
    #browser()
  }
  optA1 <- apply(EYoptMatrix,1,which.max)
  m.opt.stage1 <- colMeans(EYoptMatrix)
  return(list(m.opt.stage1=m.opt.stage1,optA1=optA1))
  
}

OptTrt.sim2b.emplik <- function(etaHMC, omegaHMC, x1, XiMatrix){
  
  B <- nrow(etaHMC)
  eta1HMC <- etaHMC[,1:(ncol(etaHMC)/2)]
  eta0HMC <- etaHMC[,(ncol(etaHMC)/2+1):(ncol(etaHMC))]
  A1.size <- ((ncol(etaHMC)-2)/4)
  eta13 <- eta1HMC[,ncol(eta1HMC)]
  eta03 <- eta0HMC[,ncol(eta0HMC)]
  EYoptMatrix <- matrix(0,nrow=B,ncol=A1.size)

  for(s in 1:A1.size){
    
    eta11.temp <- eta1HMC[,((s-1)*2+1)]
    eta01.temp <- eta0HMC[,((s-1)*2+1)]
    eta12.temp <- eta1HMC[,(s*2)]
    eta02.temp <- eta0HMC[,(s*2)]
    omega1.temp <- omegaHMC[,((s-1)*2+1)]
    omega2.temp <- omegaHMC[,(s*2)]
    #browser()
    #m.u <- (eta11.temp - eta01.temp) + (eta13 - eta03)*omega1.temp + ((eta12.temp - eta02.temp) + (eta13 - eta03)*omega2.temp)*x1
    #EYoptMatrix[,s] <- 0.5*( (eta11.temp + eta01.temp) + (eta13 + eta03)*omega1.temp + ( (eta12.temp + eta02.temp) + omega2.temp*(eta13 + eta03))*x1 + sqrt(sigmasqHMC)*sqrt(2/pi) * exp( -m.u^2/(2*sigmasqHMC) ) + m.u*( 1- 2*pnorm(-m.u/sqrt(sigmasqHMC)) ) )
    frontVec <- 0.5*((eta11.temp + eta01.temp) + (eta12.temp + eta02.temp)*x1 + (eta13 + eta03)*(omega1.temp + omega2.temp*x1))
    #browser()
    backVec <- rowMeans(apply( XiMatrix, 2, function(s){ abs(0.5*((eta13 - eta03)*( omega1.temp + omega2.temp*x1 + s ) + (eta11.temp - eta01.temp) + (eta12.temp - eta02.temp)*x1))   }  ))
    #browser()
    EYoptMatrix[,s] <- frontVec + backVec
  }
  optA1 <- apply(EYoptMatrix,1,which.max)
  m.opt.stage1 <- colMeans(EYoptMatrix)
  return(list(m.opt.stage1=m.opt.stage1,optA1=optA1))
  
}



DistributedComp.ctsX.sim2b  <- function(Y, x2, stage1design, stage2design, etaHat, omegaHat, RSS.Y, RSS.X2, YorX2 ){
  
  n.train <- length(Y)
  retObj <- NULL
  LMdim <- ncol(stage2design)
  if(YorX2==1){
    
    W <- solve(t(stage2design)%*%stage2design)
    Tsamples <- rmvt(n = 5000, sigma = ((RSS.Y/(n.train-LMdim))*W), df = (n.train-LMdim), delta = etaHat, type = "shifted")
    retObj <- list(Tsamples=Tsamples)
  }
  else{
    
    W.X2 <- solve(t(stage1design)%*%stage1design)
    LMX2dim <- ncol(stage1design)
    Tsamples.X2 <- rmvt(n = 5000, sigma = ((RSS.X2[1]/(n.train-LMX2dim))*W.X2), df = (n.train-LMX2dim), delta = omegaHat, type = "shifted")
    Res.omega.1 <- apply(Tsamples.X2,1,function(v){ sum((x2-c(stage1design %*%v))^2) })
    sigma2g.samples <- 1/rgamma(5000,shape=rep((n.train/2),5000), rate= (0.5*Res.omega.1))
    
    
    retObj <- list(Tsamples.X2=Tsamples.X2, sigma2g.samples=sigma2g.samples)
  }
  return(retObj)
  
  
}

DistributedComp.ctsX.sim2b.Emplik  <- function(y, x2, stage1design, stage2design, etaHat, omegaHat, YorX2 ,nsamp ){
  
  n.train <- length(y)
  retObj <- NULL
  LMdim <- ncol(stage2design)
  if(YorX2==1){
    
    HMC.emplik.Y <- ELHMC(initial = unname(etaHat), data = unname(cbind(stage2design,y)), FUN = FSLR, DFUN = DFSLR, n.samples = nsamp, prior = SLRpriorC, dprior = log_SLRprior_gradient_SLRpriorC, epsilon=0.005, lf.steps = 20)
    retObj <- list(accept.rate=HMC.emplik.Y$acceptance.rate, Tsamples = HMC.emplik.Y$samples)
    
  }
  else{
    
    HMC.emplik.X2 <- ELHMC(initial = unname(omegaHat), data = unname(cbind(stage1design,x2)), FUN = FSLR, DFUN = DFSLR, n.samples = nsamp, prior = SLRpriorC, dprior = log_SLRprior_gradient_SLRpriorC, epsilon=0.01, lf.steps = 20)
    retObj <- list(accept.rate=HMC.emplik.X2$acceptance.rate, Tsamples = HMC.emplik.X2$samples)
  }
  return(retObj)
  
  
}

ComputeVhat <- function(vy, X, lambda, theta){
  
  n <- length(vy)
  p <- ncol(X)
  H <- GSLRcpp(theta, cbind(X,vy))
  UnnormWeightsEL <- c(1/((1 + H %*%lambda)))
  
  SumNablaH <- matrix(0,nrow=p,ncol=p)
  for(i in 1:n){
    
    SumNablaH <- SumNablaH - matrix(X[i,],ncol=1) %*% matrix(X[i,],nrow=1)/(UnnormWeightsEL[i])
    
  }
  MiddlePart <- solve(t(H) %*% H)
  Vhat <- (t(X)%*%X) %*%MiddlePart%*%(t(X)%*%X)
  
  return(Vhat)
  
}


MCMCcNormal <- function(vy, X, B=500000, theta.init = rep(0,ncol(X))){
  
  n <- length(vy)
  p <- ncol(X)
  Xy <- cbind(X,vy)
  
  theta.MELE <- coef(lm(vy~X[,-1]))
  ELMELEobj <- evalGel(g=GSLRcpp,x=Xy,tet0=theta.MELE, optlam="nlminb")
  #browser()
  Vhat <- ComputeVhat(vy = vy, X=X, lambda = c(-ELMELEobj$lambda), theta = theta.MELE)
  invVhat <- solve(Vhat)
  thetaLibrary <- rmvnorm(n=B, mean = theta.MELE, sigma = invVhat)
  logFlipLibrary <- log(runif(B))
  theta.curr <- theta.init
  theta.new <- NULL
  thetaStore <- matrix(0,nrow=B,ncol=p)
  ELobj.curr <- evalGel(g=GSLRcpp,x=Xy,tet0=theta.curr, optlam="nlminb")
  logEL.curr <- sum(log(ELobj.curr$pt))
  logPost.curr <- sum(dnorm(theta.curr,mean=0,sd=100,log = TRUE)) + logEL.curr
  accept.count <- 0
  #browser()
  for(it in 1:B){
    
    theta.new <- thetaLibrary[it,]
    Gobjnew <- GSLRcpp(theta.new, Xy)
    lambdaObj.new <- getLamb(gt=Gobjnew,type="EL",method="optim")
    
    if(lambdaObj.new$convergence$convergence==0){
      
      lambda.curr <- -lambdaObj.new$lambda
      logEL.new <- sum(-log((1 + c(Gobjnew %*% lambda.curr)))) - n*log(n)
      
    }
    else{
      
      logEL.new <- -1000^20
      
    }
    logPost.new <- sum(dnorm(theta.new,mean=0,sd=100,log=TRUE)) + logEL.new
    logQratio <- dmvnorm(theta.curr, mean = theta.MELE, sigma = invVhat, log=TRUE) - dmvnorm(theta.new, mean = theta.MELE, sigma = invVhat, log=TRUE)
    alpha <- logPost.new - logPost.curr + logQratio
    #browser()
    logFlip <- logFlipLibrary[it]
    if(logFlip < alpha){
      
      theta.curr <- theta.new
      logEL.curr <- logEL.new
      logPost.curr <- logPost.new
      accept.count <- accept.count + 1
      
    }
    
    thetaStore[it,] <- theta.curr
    #browser()
    if(it%%100==0){
      
      cat("Completed", it, "\n")
      #gc()
      
    }
    #gc()
    
  }
  
  accept.rate <- accept.count/B
  
  return(list(accept.rate=accept.rate, thetaStore=thetaStore))
  
}

CalcAsymVar <- function(vy, X, theta.eval){
  
  n <- length(vy)
  eps <- c(vy - c(X%*%theta.eval))
  D <- - t(X)%*%X
  S <- t(X)%*%diag(eps^2)%*%X
  Vhat <- (1/n)*t(D)%*%psolve( a = S, b =  D )
  #browser()
  return((1/n)*psolve(a=Vhat, b= diag(1,ncol(Vhat)) ))
}

ChAELSVBNormal <- function(vy, X, a=(0.000001 * max(1,0.5*log(length(vy)))), mu.init=NULL, maxRuns=100001){
  
  p <- length(mu.init)
  n <- length(vy)
  Xy <- cbind(X,vy)
  Xty <- t(X)%*%vy
  XtX <- t(X)%*%X
  
  
  zMAT <- matrix( rnorm(p*maxRuns), nrow=maxRuns, ncol=p)
  mu.curr <- mu.init
  C.curr <- t(chol(CalcAsymVar(vy=vy, X = X, theta.eval = mu.curr) + diag(0.0001,length(mu.curr)) ))

  rho <- 0.9
  epsilon <- 10^(-9)
  
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mustore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim = c(p,p,maxRuns))
  
  Accugrad.mu <- rep(0,p)
  Change.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  
  LowerID <- matrix(1,nrow=p,ncol=p)
  LowerID[row(LowerID) < col(LowerID)] <- 0
  LowerID <- c(LowerID)
  
  it<-0
  maxLhat <- -9999999999999
  #browser()
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + mu.curr
    Gobj <- GSLRcpp(tet=theta.curr, x = Xy)
    Extrah <- c(-(a/n)*(Xty - XtX%*%theta.curr))
    GandExtrah <- rbind(Gobj,Extrah)
    lambda.curr <- -getLamb(gt=GandExtrah,type="EL",method="optim")$lambda
    logEL <- sum(-log(1 + c(GandExtrah %*% lambda.curr))) - (n+1)*log(n+1)
    #browser()
    grad.curr <- c(NablaAELCh2(X=X, Hmat = Gobj, lambda = lambda.curr, extrah = Extrah, a = a,tet=theta.curr))
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    grad.C <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*(grad.C^2)
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * grad.C
    C.new <- C.curr + Change.C
    C.new[row(C.new) < col(C.new)] <- 0
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    mu.curr <- mu.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mustore[it,] <- mu.curr
    
    
    #cat("Completed iter: ",it,"\n")
    
    
    
    
    #browser()
    
    
    
  }
  
  return(list(mustore = mustore, Cstore = Cstore, gradstore=gradstore))
}

ChAELSVBNormalTwoDependentVarbs <- function(vy1, vy2, X, a=(0.000001 * max(1,0.5*log(length(vy1)))), mu.init=NULL, maxRuns=100001){
  
  p <- length(mu.init)
  n <- length(vy1)
  Xy <- cbind(X,vy1,vy2)
  Xty1 <- t(X)%*%vy1
  Xty2 <- t(X)%*%vy2
  XtX <- t(X)%*%X
  
  
  zMAT <- matrix( rnorm(p*maxRuns), nrow=maxRuns, ncol=p)
  mu.curr <- mu.init
  #browser()
  C.curr1 <- t(chol(CalcAsymVar(vy=vy1, X = X, theta.eval = mu.curr[1:(p/2)]) + diag(0.0001,length(mu.curr[1:(p/2)])) ))
  C.curr2 <- t(chol(CalcAsymVar(vy=vy2, X = X, theta.eval = mu.curr[(p/2+1):p]) + diag(0.0001,length(mu.curr[(p/2+1):p])) ))
  C.curr <- blockdiag(C.curr1,C.curr2)
  
  rho <- 0.9
  epsilon <- 10^(-9)
  
  gradstore <- matrix(0,nrow=maxRuns,ncol=p)
  mustore <- matrix(0,nrow=maxRuns,ncol=p)
  Cstore <- array(0,dim = c(p,p,maxRuns))
  
  Accugrad.mu <- rep(0,p)
  Change.mu <- rep(0,p)
  Accuchange.mu <- rep(0,p)
  Accugrad.C <- matrix(0,nrow=p,ncol=p)
  Accuchange.C <- matrix(0,nrow=p,ncol=p)
  
  LowerID <- matrix(1,nrow=p,ncol=p)
  LowerID[row(LowerID) < col(LowerID)] <- 0
  LowerID <- c(LowerID)
  
  it<-0
  maxLhat <- -9999999999999
  #browser()
  
  for(it in 1:maxRuns){
    
    zdraw <- zMAT[it,]
    theta.curr <- c(C.curr %*% zdraw) + mu.curr
    #browser()
    Gobj <- STARDIntermediateRcpp(tet=theta.curr, x = Xy)
    Extrah <- c(-(a/n)*colSums(Gobj))
    GandExtrah <- rbind(Gobj,Extrah)
    lambda.curr <- -getLamb(gt=GandExtrah,type="EL",method="optim")$lambda
    logEL <- sum(-log(1 + c(GandExtrah %*% lambda.curr))) - (n+1)*log(n+1)
    #browser()
    grad.curr <- c(NablaAELChSTARD(X = X, Hmat = Gobj, lambda = lambda.curr, extrah = Extrah, a = a, tet = theta.curr))
    Accugrad.mu <- rho*Accugrad.mu + (1-rho)*grad.curr^2
    Change.mu <- sqrt(Accuchange.mu +  epsilon)/sqrt(Accugrad.mu + epsilon) *grad.curr
    mu.new <- mu.curr + Change.mu
    Accuchange.mu <- rho*Accuchange.mu + (1-rho)*Change.mu^2
    
    grad.C <- (matrix(grad.curr,ncol=1) %*% matrix(zdraw,nrow=1) + diag(1/diag(C.curr)) )
    Accugrad.C <- rho*Accugrad.C  +(1-rho)*(grad.C^2)
    Change.C <- sqrt(Accuchange.C +  epsilon )/ sqrt( Accugrad.C + epsilon ) * grad.C
    C.new <- C.curr + Change.C
    C.new[row(C.new) < col(C.new)] <- 0
    Accuchange.C <- rho*Accuchange.C + (1-rho)*Change.C^2
    Accuchange.C[row(Accuchange.C)< col(Accuchange.C)] <- 0
    
    mu.curr <- mu.new
    C.curr <- C.new
    gradstore[it,] <- grad.curr
    Cstore[,,it] <- C.curr
    mustore[it,] <- mu.curr
    
    # if((it%%10)==0){
    #   
    #   cat("Completed iter: ",it,"\n")
    #   
    # }
    #
    
    
    
    
    #browser()
    
    
    
  }
  
  return(list(mustore = mustore, Cstore = Cstore, gradstore=gradstore))
}


DistributedComp.ctsX.sim2b.Emplik.VB  <- function(y, x2, stage1design, stage2design, etaHat, omegaHat, YorX2 ,nsamp ){
  
  n.train <- length(y)
  retObj <- NULL
  LMdim <- ncol(stage2design)
  if(YorX2==1){
    
    #HMC.emplik.Y <- ELHMC(initial = unname(etaHat), data = unname(cbind(stage2design,y)), FUN = FSLR, DFUN = DFSLR, n.samples = nsamp, prior = SLRpriorC, dprior = log_SLRprior_gradient_SLRpriorC, epsilon=0.005, lf.steps = 20)
    VB.emplik.Y <- ChAELSVBNormal(vy = y, X = stage2design, mu.init = etaHat, maxRuns = 10000)
    retObj <- list(Tsamples = rmvnorm(nsamp, mean = VB.emplik.Y$mustore[10000,], sigma = (VB.emplik.Y$Cstore[,,10000]%*%t(VB.emplik.Y$Cstore[,,10000])) ) )
    
  }
  else{
    
    #HMC.emplik.X2 <- ELHMC(initial = unname(omegaHat), data = unname(cbind(stage1design,x2)), FUN = FSLR, DFUN = DFSLR, n.samples = nsamp, prior = SLRpriorC, dprior = log_SLRprior_gradient_SLRpriorC, epsilon=0.01, lf.steps = 20)
    VB.emplik.X2 <- ChAELSVBNormal(vy = x2, X = stage1design, mu.init = omegaHat, maxRuns = 10000)
    retObj <- list(Tsamples = rmvnorm(nsamp, mean = VB.emplik.X2$mustore[10000,], sigma = (VB.emplik.X2$Cstore[,,10000]%*%t(VB.emplik.X2$Cstore[,,10000])) ))
  }
  return(retObj)
  
  
}


ql.glm.ctsX.sim2b = function(dataset, o1.test, stage2design, stage1design){
  # data
  n = nrow(dataset$data); o1 = dataset$data[,"o1"]; a1 = dataset$data[,"a1"]
  o2 = dataset$data[,"o2"]; a2 = dataset$data[,"a2"]; y = dataset$data[,"y"]
  n.test <- length(o1.test)
  
  # Stage 2 Parameter Estimation and Percentile Bootstrap
  X2tX2.inv = solve(t(stage2design)%*%stage2design); beta2.hat = X2tX2.inv%*%t(stage2design)%*%y
  mu2.hat = as.vector(stage2design%*%beta2.hat)
  stage2design.counterfact <- cbind(stage2design[,(ncol(stage2design)/2 + 1):ncol(stage2design)], stage2design[,1:(ncol(stage2design)/2)])
  y.tilde <- pmax(mu2.hat, c(stage2design.counterfact%*%beta2.hat))
  
  
  # Stage 1 Parameter Estimation
  X1tX1.inv = solve(t(stage1design)%*%stage1design); beta1.hat = c(X1tX1.inv%*%t(stage1design)%*%y.tilde)
  stage1pred <- cbind(1,o1.test)
  A1.size <- ncol(stage1design)/2
  mu1.hat <- NULL
  for(k in 1:A1.size){
    
    mu1.hat <- cbind(mu1.hat, c(stage1pred%*% beta1.hat[ ((k-1)*2+1):(2*k)]))
    
  }
  a1opt <- apply(mu1.hat,1,which.max)
  
  return(list(beta2.hat=beta2.hat, beta1.hat=beta1.hat, a1opt=a1opt,mu1.hat=mu1.hat))
}

FSLR <- function(params, X){
  
  p <- ncol(X) - 1
  
  beta <- params
  y <- X[,(p+1)]
  x.predictors <- X[,1:p]
  u <- y - c(x.predictors %*% beta)
  ans <- matrix(rep(u,p)*c(x.predictors),nrow=nrow(X),ncol=p)
  return(ans)
  
}

DFSLR <- function(params, X){
  
  p <- ncol(X) - 1
  n <- nrow(X)
  beta <- params
  y <- X[,(p+1)]
  x.predictors <- X[,1:p]
  #u <- y - c(x.predictors %*% beta)
  CrossProd <- apply(x.predictors,1,function(s){ matrix(s,ncol=1) %*% matrix(s,nrow=1) })
  #ansvec <- c(CrossProd) * rep(u,rep(p*p,n))
  ansvec <- -c(CrossProd)
  ans <- array(ansvec,dim = c(p,p,n))
  return(ans)
  
}

SLRpriorC <- function(x){
  
  p <- length(x)
  kappa <- 100
  exp(sum( dnorm(x,mean=0, sd = kappa,log = TRUE) ))
  
  
  
}

log_SLRprior_gradient_SLRpriorC <- function(x){
  
  kappa <- 100
  -x/(kappa^2)
  
  
}

FSLRdouble <- function(params, X){
  
  p <- ncol(X) - 2
  
  beta1 <- params[1:p]
  beta2 <- params[(p+1):(2*p)]
  y1 <- X[,(p+1)]
  y2 <- X[,(p+2)]
  x.predictors <- X[,1:p]
  u1 <- y1 - c(x.predictors %*% beta1)
  u2 <- y2 - c(x.predictors %*% beta2)
  ans <- cbind( matrix(rep(u1,p)*c(x.predictors),nrow=nrow(X),ncol=p), matrix( rep(u2,p)*c(x.predictors), nrow=nrow(X), ncol=p) )
  return(ans)
  
}

DFSLRdouble <- function(params, X){
  
  p <- ncol(X) - 2
  n <- nrow(X)
  
  beta1 <- params[1:p]
  beta2 <- params[(p+1):(2*p)]
  x.predictors <- X[,1:p]
  
  CrossProd <- apply(x.predictors,1,function(s){ MS <- matrix(s,ncol=1) %*% matrix(s,nrow=1); as.matrix(bdiag(MS,MS)) })
  #ansvec <- c(CrossProd) * rep(u,rep(p*p,n))
  ansvec <- -c(CrossProd)
  ans <- array(ansvec,dim = c(2*p,2*p,n))
  return(ans)
  
}

DistributedComp.ctsX.sim3  <- function(Y, x2, stage1design, stage2design, etaHat, omegaHat, RSS.Y, RSS.X2, YorX2 ){
  
  n.train <- length(Y)
  retObj <- NULL
  LMdim <- ncol(stage2design)
  if(YorX2==1){
    
    W <- psolve(t(stage2design)%*%stage2design)
    Tsamples <- rmvt(n = 5000, sigma = ((RSS.Y/(n.train-LMdim))*W), df = (n.train-LMdim), delta = etaHat, type = "shifted")
    retObj <- list(Tsamples=Tsamples)
  }
  else{
    
    W.X2 <- psolve(t(stage1design)%*%stage1design)
    LMX2dim <- ncol(stage1design)
    Tsamples.X2 <- rmvt(n = 5000, sigma = ((RSS.X2[1]/(n.train-LMX2dim))*W.X2), df = (n.train-LMX2dim), delta = omegaHat, type = "shifted")
    Res.omega.1 <- apply(Tsamples.X2,1,function(v){ sum((x2-c(stage1design %*%v))^2) })
    sigma2g.samples <- 1/rgamma(5000,shape=rep((n.train/2),5000), rate= (0.5*Res.omega.1))
    
    
    retObj <- list(Tsamples.X2=Tsamples.X2, sigma2g.samples=sigma2g.samples)
  }
  return(retObj)
  
  
}

OptTrt.ctsX.sim3 <- function(eta1MC, eta0MC, omegaMC, sigma2gMC, KnotsMatrix, x1.test.point, x1.basis){
  
  B <- nrow(eta1MC)
  Bprime <- 50
  repx1test <- rep(x1.test.point,B)
  Emu.draws <- matrix(0,nrow=2,ncol=B)
  
  for(s in 1:2){
    
    eta1MC.reduced <- eta1MC[,c(((s-1)*ncol(eta1MC)/2+1):((s*ncol(eta1MC))/2))]
    eta0MC.reduced <- eta0MC[,c(((s-1)*ncol(eta1MC)/2+1):((s*ncol(eta1MC))/2))]
    omegaMC.reduced <-  omegaMC[,((s-1)*ncol(omegaMC)/2+1):( (s*ncol(omegaMC))/2) ]
    max.draws <- matrix(0,nrow=Bprime,ncol=B)
    common.mean <- c(omegaMC.reduced%*%x1.basis)
    X2draws.BIG <- rep(common.mean,Bprime) + rnorm( (B*Bprime) ,sd=rep(sqrt(sigma2gMC),Bprime))
    X2BASIS.BIG.temp <- basis.tps( x=cbind(rep(repx1test,Bprime),X2draws.BIG),  knots = KnotsMatrix, intercept = TRUE  )
    max.draws <- matrix(0,nrow=Bprime,ncol=B)
    for(temp.id in 1:Bprime){

      x2.basis <- X2BASIS.BIG.temp[((temp.id -1)*B+1):(temp.id*B),]
      max.draws[temp.id,] <- pmax(apply(eta1MC.reduced*x2.basis,1,sum),apply(eta0MC.reduced*x2.basis,1,sum))
      
    }
    Emu.draws[s,] <- colMeans(max.draws)
    
  }
  
  optA1 <- c(apply(Emu.draws,2,function(r){ which.max(r) }))
  Emu.means <- rowMeans(Emu.draws)
  
  return(list(optA1=optA1,Emu.means=Emu.means) )
  
}

OptTrt.sim3.emplik <- function(eta1MC, eta0MC, omegaMC, KnotsMatrix, x1.test.point, x1.basis, XiMatrix){
  
  B <- nrow(eta1MC)
  repx1test <- rep(x1.test.point,B)
  Emu.draws <- matrix(0,nrow=3,ncol=B)
  nTrain <- ncol(XiMatrix)
  Q1est <- array(0,dim = c(B,nTrain,2))
  #browser()
  for(s in 1:2){
    
    eta1MC.reduced <- eta1MC[,c(((s-1)*ncol(eta1MC)/2+1):((s*ncol(eta1MC))/2))]
    #browser()
    eta0MC.reduced <- eta0MC[,c(((s-1)*ncol(eta1MC)/2+1):((s*ncol(eta1MC))/2))]
    #browser()
    omegaMC.reduced <-  omegaMC[,((s-1)*ncol(omegaMC)/2+1):( (s*ncol(omegaMC))/2) ]
    #browser()
    common.mean <- c(omegaMC.reduced%*%x1.basis)
    #browser()
    X2.discretized <- apply(XiMatrix,2,function(Xi){ Xi + common.mean })
    #Use browser to check dimension of aboe output and then flatten into a B*Bprime vector to sum with rep(common.mean,Bprime)
    #X2BASIS.BIG.temp <- basis.tps(cbind(x1.test.point,c(t(X2.discretized)) ), knots = KnotsMatrix, intercept = TRUE  )
    #browser()
    X2BASIS.BIG.temp <- basis.tps(cbind(x1.test.point,c(X2.discretized) ), knots = KnotsMatrix, intercept = TRUE  )
    for(temp.id in 1:nTrain){
      
      #browser()
      Q1est[,temp.id,s] <- pmax( rowSums(X2BASIS.BIG.temp[((temp.id-1)*B +1):(temp.id*B),]*eta1MC.reduced) ,  rowSums(X2BASIS.BIG.temp[((temp.id-1)*B +1):(temp.id*B),]*eta0MC.reduced) )
      #browser()
    }
    rm(X2BASIS.BIG.temp)
    gc();gc();gc()
    
  }
  optA1 <- apply(rbind(rowMeans(Q1est[,,1]),rowMeans(Q1est[,,2])),2,which.max)
  m.opt.stage1 <- c(mean(Q1est[,,1]),mean(Q1est[,,2]))
  return(list(m.opt.stage1=m.opt.stage1,optA1=optA1))
  
}

ql.sim3 = function(dataset, o1.test, stage1.test.basis){
  # data
  n = nrow(dataset$data); o1 = dataset$data[,"o1"]; a1 = dataset$data[,"a1"]
  o2 = dataset$data[,"o2"]; a2 = dataset$data[,"a2"]; y = dataset$data[,"y"]
  stage1.basis = dataset$stage1basis; stage2.basis = dataset$stage2basis
  n.test <- length(o1.test)
  #browser()
  # implement method
  stage2.pred <- cbind(apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==1)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==1)}), apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==0)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==0)}) )
  #browser()
  beta2.hat <- psolve(t(stage2.pred)%*%stage2.pred, t(stage2.pred)%*%y)
  stage2.pred.counterfact <- cbind(apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==0)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==0)}), apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==1)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==1)}) )
  #browser()
  ytilde <- pmax(c(stage2.pred%*%beta2.hat), c(stage2.pred.counterfact%*%beta2.hat))
  stage1.pred <- cbind(apply(stage1.basis,2,function(s){ s*(a1==-1)}), apply(stage1.basis,2,function(s){ s*(a1==1)}) )
  beta1.hat <- c(psolve(t(stage1.pred)%*%stage1.pred, t(stage1.pred)%*%ytilde))
  #browser()
  EYopt.by.A1 <- rbind( c(stage1.test.basis %*% beta1.hat[1:(ncol(stage1.pred)/2)]), c(stage1.test.basis %*% beta1.hat[(ncol(stage1.pred)/2 + 1):(ncol(stage1.pred))])  )
  a1.opt <- apply(EYopt.by.A1,2,which.max)
  mu1.hat <- apply(EYopt.by.A1,2,max)
  #browser()
  
  return(list(beta2.hat=beta2.hat,beta1.hat=beta1.hat,ytilde=ytilde,mu1.hat=mu1.hat, a1.opt=a1.opt, EYopt.by.A1=EYopt.by.A1))
}




ql.sim4 = function(dataset, o1.test, stage1.test.basis){
  
  n = nrow(dataset$data); o1 = dataset$data[,"o1"]; a1 = dataset$data[,"a1"]
  o2 = dataset$data[,"o2"]; a2 = dataset$data[,"a2"]; y = dataset$data[,"y"]
  stage1.basis = dataset$stage1basis; stage2.basis = dataset$stage2basis
  n.test <- length(o1.test)
  
  stage2.pred <- cbind(apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==1)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==1)}), apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==0)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==0)}))
  beta2.hat <- psolve(t(stage2.pred)%*%stage2.pred, t(stage2.pred)%*%y)
  stage2.pred.counterfact <- cbind(apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==0)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==0)}), apply(stage2.basis,2,function(s){ s*(a1==-1)*(a2==1)}), apply(stage2.basis,2,function(s){ s*(a1==1)*(a2==1)}))
  ytilde <- pmax(c(stage2.pred%*%beta2.hat), c(stage2.pred.counterfact%*%beta2.hat))
  stage1.pred <- cbind(apply(stage1.basis,2,function(s){ s*(a1==-1)}), apply(stage1.basis,2,function(s){ s*(a1==1)}) )
  beta1.hat <- c(psolve(t(stage1.pred)%*%stage1.pred, t(stage1.pred)%*%ytilde))
  mu1.hat <- c(c(stage1.test.basis %*% beta1.hat[1:(ncol(stage1.pred)/2)]),c(stage1.test.basis %*% beta1.hat[(ncol(stage1.pred)/2 + 1):ncol(stage1.pred)]))
  a1.opt <- sign(mu1.hat[(n.test+1):(2*n.test)] - mu1.hat[1:n.test])
  
  
  return(list(beta2.hat=beta2.hat,beta1.hat=beta1.hat,ytilde=ytilde,mu1.hat=mu1.hat, a1.opt=a1.opt))
}

DistributedComp.ctsX.sim4  <- function(Y, x2, stage1design, stage2design, etaHat, omegaHat, RSS.Y, RSS.X2, YorX2 ){
  
  n.train <- length(Y)
  retObj <- NULL
  LMdim <- ncol(stage2design)
  if(YorX2==1){
    
    W <- psolve(t(stage2design)%*%stage2design)
    Tsamples <- rmvt(n = 5000, sigma = ((RSS.Y/(n.train-LMdim))*W), df = (n.train-LMdim), delta = etaHat, type = "shifted")
    retObj <- list(Tsamples=Tsamples)
  }
  else{
    
    W.X2 <- psolve(t(stage1design)%*%stage1design)
    LMX2dim <- ncol(stage1design)
    Tsamples.X2 <- rmvt(n = 5000, sigma = ((RSS.X2[1]/(n.train-LMX2dim))*W.X2), df = (n.train-LMX2dim), delta = omegaHat, type = "shifted")
    Res.omega.1 <- apply(Tsamples.X2,1,function(v){ sum((x2-c(stage1design %*%v))^2) })
    sigma2g.samples <- 1/rgamma(5000,shape=rep((n.train/2),5000), rate= (0.5*Res.omega.1))
    
    
    retObj <- list(Tsamples.X2=Tsamples.X2, sigma2g.samples=sigma2g.samples)
  }
  return(retObj)
  
  
}

OptTrt.ctsX.sim4 <- function(eta1MC, eta0MC, omegaMC, sigma2gMC, KnotsMatrix, x1.test.point, x1.basis){
  
  B <- nrow(eta1MC)
  Bprime <- 50
  repx1test <- rep(x1.test.point,B)
  ######Stage 1 when s=1######
  eta1MC.reduced <- eta1MC[,c((ncol(eta1MC)/2+1):ncol(eta1MC))]
  eta0MC.reduced <- eta0MC[,c((ncol(eta1MC)/2+1):ncol(eta1MC))]
  omegaMC.reduced <-  omegaMC[,c((ncol(omegaMC)/2+1):ncol(omegaMC))]
  common.mean <- c(omegaMC.reduced%*%x1.basis)
  max.draws <- matrix(0,nrow=Bprime,ncol=B)
  X2draws.BIG <- rep(common.mean,Bprime) + rnorm( (B*Bprime) ,sd=rep(sqrt(sigma2gMC),Bprime))
  X2BASIS.BIG.temp <- basis.tps( x=cbind(rep(repx1test,Bprime),X2draws.BIG),  knots = KnotsMatrix, intercept = TRUE  )
  for(temp.id in 1:Bprime){
    
    x2.basis <- X2BASIS.BIG.temp[((temp.id -1)*B+1):(temp.id*B),]
    max.draws[temp.id,] <- pmax(apply(eta1MC.reduced*x2.basis,1,sum),apply(eta0MC.reduced*x2.basis,1,sum))
    
  }
  Emubar.s1.draws <- colMeans(max.draws)
  
  ######Stage 1 when s=-1######
  eta1MC.reduced <- eta1MC[,c(1:(ncol(eta1MC)/2))]
  eta0MC.reduced <- eta0MC[,c(1:(ncol(eta1MC)/2))]
  omegaMC.reduced <-  omegaMC[,c(1:(ncol(omegaMC)/2))]
  common.mean <- c(omegaMC.reduced%*%x1.basis)
  X2draws.BIG <- rep(common.mean,Bprime) + rnorm( (B*Bprime) ,sd=rep(sqrt(sigma2gMC),Bprime))
  X2BASIS.BIG.temp <- basis.tps( x=cbind(rep(repx1test,Bprime),X2draws.BIG),  knots = KnotsMatrix, intercept = TRUE  )
  max.draws <- matrix(0,nrow=Bprime,ncol=B)
  
  for(temp.id in 1:Bprime){
    
    x2.basis <- X2BASIS.BIG.temp[((temp.id -1)*B+1):(temp.id*B),]
    max.draws[temp.id,] <- pmax(apply(eta1MC.reduced*x2.basis,1,sum),apply(eta0MC.reduced*x2.basis,1,sum))
    
  }
  Emubar.sNeg1.draws <- colMeans(max.draws)
  
  optA1.vec <- sign(Emubar.s1.draws - Emubar.sNeg1.draws)
  optA1 <- getmode(optA1.vec)
  
  #Output dimensions of list items: number of MC samples (B) by number of test data (nTest)
  return(list(optA1=optA1,m1opt.s1=mean(Emubar.s1.draws),m1opt.sNeg1=mean(Emubar.sNeg1.draws) ))
  
}

OptTrt.sim4.emplik <- function(eta1MC, eta0MC, omegaMC, KnotsMatrix, x1.test.point, x1.basis, XiMatrix){
  
  B <- nrow(eta1MC)
  repx1test <- rep(x1.test.point,B)
  Emu.draws <- matrix(0,nrow=3,ncol=B)
  nTrain <- ncol(XiMatrix)
  Q1est <- array(0,dim = c(B,nTrain,2))
  m.opt.stage1 <- matrix(0,nrow=2,ncol=B)

  for(s in 1:2){
    
    eta1MC.reduced <- eta1MC[,c(((s-1)*ncol(eta1MC)/2+1):((s*ncol(eta1MC))/2))]
    eta0MC.reduced <- eta0MC[,c(((s-1)*ncol(eta1MC)/2+1):((s*ncol(eta1MC))/2))]
    omegaMC.reduced <-  omegaMC[,((s-1)*ncol(omegaMC)/2+1):( (s*ncol(omegaMC))/2) ]
    common.mean <- c(omegaMC.reduced%*%x1.basis)
    X2.discretized <- apply(XiMatrix,2,function(Xi){ Xi + common.mean })
    #X2BASIS.BIG.temp <- basis.tps(cbind(x1.test.point,c(X2.discretized) ), knots = KnotsMatrix, intercept = TRUE  )
    optA1 <- rep(0,nTrain)
    for(temp.id in 1:nTrain){
      
      X2BASIS.BIG.temp <- basis.tps(cbind(x1.test.point,X2.discretized[,temp.id] ), knots = KnotsMatrix, intercept = TRUE  )
      #Q1est[,temp.id,s] <- pmax( rowSums(X2BASIS.BIG.temp[((temp.id-1)*B +1):(temp.id*B),]*eta1MC.reduced) ,  rowSums(X2BASIS.BIG.temp[((temp.id-1)*B +1):(temp.id*B),]*eta0MC.reduced) )
      Q1est[,temp.id,s] <- pmax( rowSums(X2BASIS.BIG.temp*eta1MC.reduced) ,  rowSums(X2BASIS.BIG.temp*eta0MC.reduced) )
      #pmax(rowSums(X2BASIS.BIG.temp*eta1MC.reduced),rowSums(X2BASIS.BIG.temp*eta0MC.reduced))
      
    }
    rm(X2BASIS.BIG.temp)
    gc();gc();gc()
    
  }
  optA1 <- apply(rbind(rowMeans(Q1est[,,1]),rowMeans(Q1est[,,2])),2,which.max)
  m.opt.stage1 <- c(mean(Q1est[,,1]),mean(Q1est[,,2]))
  return(list(m.opt.stage1=m.opt.stage1,optA1=optA1))
  
}

DistributedComp.ctsX.sim3.Emplik.VB  <- function(y, x2, stage1design, stage2design, etaHat, omegaHat, YorX2 ,nsamp, VBiter ){
  
  n.train <- length(y)
  retObj <- NULL
  LMdim <- ncol(stage2design)
  if(YorX2==1){
    
    VB.emplik.Y <- ChAELSVBNormal(vy = y, X = stage2design, mu.init = etaHat, maxRuns = VBiter)
    retObj <- list(Tsamples = rmvnorm(nsamp, mean = VB.emplik.Y$mustore[VBiter,], sigma = (VB.emplik.Y$Cstore[,,VBiter]%*%t(VB.emplik.Y$Cstore[,,VBiter])) ) )
    
  }
  else{
    
    VB.emplik.X2 <- ChAELSVBNormal(vy = x2, X = stage1design, mu.init = omegaHat, maxRuns = VBiter)
    retObj <- list(Tsamples = rmvnorm(nsamp, mean = VB.emplik.X2$mustore[VBiter,], sigma = (VB.emplik.X2$Cstore[,,VBiter]%*%t(VB.emplik.X2$Cstore[,,VBiter])) ))
  }
  return(retObj)
  
  
}

DistributedComp.ctsX.sim4.Emplik.VB  <- function(y, x2, stage1design, stage2design, etaHat, omegaHat, YorX2 ,nsamp, VBiter ){
  
  n.train <- length(y)
  retObj <- NULL
  LMdim <- ncol(stage2design)
  if(YorX2==1){
    
    VB.emplik.Y <- ChAELSVBNormal(vy = y, X = stage2design, mu.init = etaHat, maxRuns = VBiter)
    retObj <- list(Tsamples = rmvnorm(nsamp, mean = VB.emplik.Y$mustore[VBiter,], sigma = (VB.emplik.Y$Cstore[,,VBiter]%*%t(VB.emplik.Y$Cstore[,,VBiter])) ) )
    
  }
  else{
    
    VB.emplik.X2 <- ChAELSVBNormal(vy = x2, X = stage1design, mu.init = omegaHat, maxRuns = VBiter)
    retObj <- list(Tsamples = rmvnorm(nsamp, mean = VB.emplik.X2$mustore[VBiter,], sigma = (VB.emplik.X2$Cstore[,,VBiter]%*%t(VB.emplik.X2$Cstore[,,VBiter])) ))
  }
  return(retObj)
  
  
}


# OptTrt.STARD.emplik <- function(eta1HMC, eta0HMC, omega1HMC, omega2HMC, x1, XiArray){
#   
#   B <- nrow(etaHMC)
#   A1.size <- 2
#   eta13 <- eta1HMC[,((ncol(eta1HMC)-1):ncol(eta1HMC))]
#   eta03 <- eta0HMC[,((ncol(eta0HMC)-1):ncol(eta0HMC))]
#   EYoptMatrix <- matrix(0,nrow=B,ncol=A1.size)
#   
#   #when s=1
#     eta11.temp <- eta1HMC[,1:3]
#     eta01.temp <- eta0HMC[,4:6]
#     eta12.temp <- eta1HMC[,4:5]
#     eta02.temp <- eta0HMC[,7:8]
#     omega11.temp <- omega1HMC[,4]
#     omega12.temp <- omega1HMC[,5:6]
#     omega21.temp <- omega2HMC[,4]
#     omega22.temp <- omega2HMC[,5:6]
#     #browser()
#     #m.u <- (eta11.temp - eta01.temp) + (eta13 - eta03)*omega1.temp + ((eta12.temp - eta02.temp) + (eta13 - eta03)*omega2.temp)*x1
#     #EYoptMatrix[,s] <- 0.5*( (eta11.temp + eta01.temp) + (eta13 + eta03)*omega1.temp + ( (eta12.temp + eta02.temp) + omega2.temp*(eta13 + eta03))*x1 + sqrt(sigmasqHMC)*sqrt(2/pi) * exp( -m.u^2/(2*sigmasqHMC) ) + m.u*( 1- 2*pnorm(-m.u/sqrt(sigmasqHMC)) ) )
#     frontVec <- 0.5*((eta11.temp + eta01.temp) + sum((eta12.temp + eta02.temp)*x1) + sum((eta13 + eta03)*c((omega11.temp + sum(omega12.temp*x1)),(omega21.temp + sum(omega22.temp*x1))))  )
#     browser()
#     backVec <- rowMeans(apply( XiArray, c(1,2), function(s){ abs(0.5*( sum((eta13 - eta03)*(c( omega11.temp + sum(omega12.temp*x1), omega21.temp + sum(omega22.temp*x1)  ) + s)  ) + (eta11.temp - eta01.temp) + sum((eta12.temp - eta02.temp)*x1)  ))   }  ))
#     #backVec <- rowMeans(apply( XiArray, 2, function(s){ abs(0.5*( sum((eta13 - eta03)*( omega1.temp + omega2.temp*x1 + s )) + (eta11.temp - eta01.temp) + sum((eta12.temp - eta02.temp)*x1)  ))   }  ))
#     browser()
#     EYoptMatrix[,s] <- frontVec + backVec
# 
#   optA1 <- apply(EYoptMatrix,1,which.max)
#   m.opt.stage1 <- colMeans(EYoptMatrix)
#   return(list(m.opt.stage1=m.opt.stage1,optA1=optA1))
#   
# }

OptTrt.STARD.emplik <- function(eta1MC, eta0MC, omega1MC, omega2MC, x1, XiMatrix1, XiMatrix2){
  
  B <- nrow(eta1MC)
  N <- ncol(XiMatrix1)
  EYoptMatrix <- matrix(0,nrow=B,ncol=2)
  ######Stage 1 when s=1######
  eta1MC.reduced <- eta1MC
  eta0MC.reduced <- eta0MC[,c(4:8)]
  omega1MC.reduced <-  omega1MC[,c(4:6)]
  omega2MC.reduced <-  omega2MC[,c(4:6)]
  lambda0 <- eta1MC.reduced[,1] + eta0MC.reduced[,1] + c((eta1MC.reduced[,2:3] + eta0MC.reduced[,2:3])%*%x1)
  #browser()
  lambda1 <- eta1MC.reduced[,1] - eta0MC.reduced[,1] + c((eta1MC.reduced[,2:3] - eta0MC.reduced[,2:3])%*%x1)
  FrontAll <- 0.5*(lambda0+rowSums((eta1MC.reduced[,4:5] + eta0MC.reduced[,4:5])*cbind(c(omega1MC.reduced[,1] + omega1MC.reduced[,2:3]%*%x1), c(omega2MC.reduced[,1] + omega2MC.reduced[,2:3]%*%x1))))
  accumu <- rep(0,B)
  #browser()
  for(i in 1:N){
    
    accumu <- accumu + 0.5*abs(lambda1 + rowSums((eta1MC.reduced[,4:5] - eta0MC.reduced[,4:5])*cbind(( omega1MC.reduced[,1] + c(omega1MC.reduced[,2:3]%*%x1) + XiMatrix1[,i]), ( omega2MC.reduced[,1] + c(omega2MC.reduced[,2:3]%*%x1) + XiMatrix2[,i]) )))
    
  }
  #browser()
  EYoptMatrix[,2] <- FrontAll + (accumu)/(N)
  
  ######Stage 1 when s=-1######
  eta0MC.reduced <- eta0MC[,c(1:3,7:8)]
  omega1MC.reduced <-  omega1MC[,c(1:3)]
  omega2MC.reduced <-  omega2MC[,c(1:3)]
  EYoptMatrix[,1] <- eta0MC.reduced[,1] + c(eta0MC.reduced[,2:3]%*%x1) + eta0MC.reduced[,4]*(omega1MC.reduced[,1] + c(omega1MC.reduced[,2:3]%*%x1)  ) + eta0MC.reduced[,5]*(omega2MC.reduced[,1] + c(omega2MC.reduced[,2:3]%*%x1)  )
  
  optA1 <- getmode(apply(EYoptMatrix,1,which.max))
  optA1 <- c(-1,1)[optA1]
  #browser()
  return(list(optA1=optA1))
  
}

STARD.DrawCoefs.Emplik.VB  <- function(y, stagedesign, estHat, nsamp ){
  
  n.train <- length(y)
  LMdim <- ncol(stagedesign)
  VB.emplik <- ChAELSVBNormal(vy = y, X = stagedesign, mu.init = estHat, maxRuns = 10000)
  retObj <- list(Tsamples = rmvnorm(nsamp, mean = VB.emplik$mustore[10000,], sigma = (VB.emplik$Cstore[,,10000]%*%t(VB.emplik$Cstore[,,10000])) ) )

  return(retObj)
}

STARD.DrawCoefs.X2.Emplik.VB  <- function(x21, x22, stage1design, omegaHat1, omegaHat2, nsamp ){
  
  n.train <- length(x21)
  omegaVB.emplik <- ChAELSVBNormalTwoDependentVarbs(vy1 = x21, vy2 = x22, X = stage1design, mu.init = c(omegaHat1,omegaHat2), maxRuns = 100 )
  retObj <- list(Tsamples = rmvnorm(nsamp, mean = omegaVB.emplik$mustore[100,], sigma = (omegaVB.emplik$Cstore[,,100]%*%t(omegaVB.emplik$Cstore[,,100])) ) )
  
  return(retObj)
}

STARD.DrawCoefs.Emplik.HMC  <- function(y, stagedesign, estHat, nsamp ){
  
  n.train <- length(y)
  LMdim <- ncol(stagedesign)
  HMC.emplik <- ELHMC(initial = unname(estHat), data = unname(cbind(stagedesign,y)), FUN = FSLR, DFUN = DFSLR, n.samples = (1000+nsamp), prior = SLRpriorC, detailed=TRUE, dprior = log_SLRprior_gradient_SLRpriorC, epsilon=0.02)
  retObj <- HMC.emplik$samples[1001:(1000+nsamp), ]
  
  return(retObj)
}

STARD.DrawCoefs.X2.Emplik.HMC  <- function(x21, x22, stage1design, omegaHat1, omegaHat2, nsamp ){
  
  n.train <- length(x21)
  omegaVB.emplik <- ELHMC(initial = unname(c(omegaHat1,omegaHat2)), data = cbind(stage1design,x21,x22), FUN = FSLR, DFUN = DFSLRdouble, n.samples = (1000+nsamp), prior = SLRpriorC, detailed=TRUE, dprior = log_SLRprior_gradient_SLRpriorC, epsilon=0.02)
  
  retObj <- NULL
  return(retObj)
}

FSLR <- function(params, X){
  
  p <- ncol(X) - 1
  
  beta <- params
  y <- X[,(p+1)]
  x.predictors <- X[,1:p]
  u <- y - c(x.predictors %*% beta)
  ans <- matrix(rep(u,p)*c(x.predictors),nrow=nrow(X),ncol=p)
  return(ans)
  
}


DFSLR <- function(params, X){
  
  p <- ncol(X) - 1
  
  beta <- params
  y <- X[,(p+1)]
  x.predictors <- X[,1:p]
  #u <- y - c(x.predictors %*% beta)
  CrossProd <- apply(x.predictors,1,function(s){ matrix(s,ncol=1) %*% matrix(s,nrow=1) })
  #ansvec <- c(CrossProd) * rep(u,rep(p*p,n))
  ansvec <- -c(CrossProd)
  ans <- array(ansvec,dim = c(p,p,n))
  return(ans)
  
}

SLRpriorC <- function(x){
  
  p <- length(x)
  kappa <- 100
  exp(sum( dnorm(x,mean=0, sd = kappa,log = TRUE) ))
  
  
  
}

log_SLRprior_gradient_SLRpriorC <- function(x){
  
  kappa <- 100
  -x/(kappa^2)
  
  
}


DistributedComp.ctsX.STARD  <- function(Y, Z2, Z1, X2, etaHat, omegaHat, RSS.Y, RSS.X2, LMdim, YorX2 , nsamps=10000){
  
  n.train <- length(Y)
  retObj <- NULL
  if(YorX2==1){
    
    W <- solve(t(Z2)%*%Z2)
    Tsamples <- rmvt(n = nsamps, sigma = ((RSS.Y/(n.train-LMdim))*W), df = (n.train-LMdim), delta = etaHat, type = "shifted")
    retObj <- list(Tsamples=Tsamples)
  }
  else{
    
    W.X2 <- solve(t(Z1)%*%Z1)
    Tsamples.X2 <- rmvt(n = nsamps, sigma = ((RSS.X2/(n.train-ncol(Z1)))*W.X2), df = (n.train-ncol(Z1)), delta = omegaHat, type = "shifted")
    Res.omega <- apply(Tsamples.X2,1,function(v){ sum((X2-c(Z1 %*%v))^2) })
    sigma2g.samples <- 1/rgamma(nsamps,shape=rep((n.train/2),nsamps), rate= (0.5*Res.omega))
    
    retObj <- list(Tsamples.X2=Tsamples.X2, sigma2g.samples=sigma2g.samples)
  }
  return(retObj)
  
  
}

OptTrt.ctsX.STARDv3 <- function(eta1MC, eta0MC, omega1MC, omega2MC, sigma2g1MC, sigma2g2MC, x1){
  
  B <- nrow(eta1MC)
  
  ######Stage 1 when s=1######
  eta1MC.reduced <- eta1MC
  eta0MC.reduced <- eta0MC[,c(4:8)]
  omega1MC.reduced <-  omega1MC[,c(4:6)]
  omega2MC.reduced <-  omega2MC[,c(4:6)]
  lambda0 <- eta1MC.reduced[,1] + eta0MC.reduced[,1] + c((eta1MC.reduced[,2:3] + eta0MC.reduced[,2:3])%*%x1)
  lambda1 <- eta1MC.reduced[,1] - eta0MC.reduced[,1] + c((eta1MC.reduced[,2:3] - eta0MC.reduced[,2:3])%*%x1)
  X2coef.diff <- cbind(eta1MC.reduced[,4] - eta0MC.reduced[,4], eta1MC.reduced[,5] - eta0MC.reduced[,5])
  X2coef.diff.sigma <- cbind(X2coef.diff[,1]*sigma2g1MC,X2coef.diff[,2]*sigma2g2MC)
  
  s.x2.sq <- X2coef.diff.sigma[,1]*X2coef.diff[,1] + X2coef.diff.sigma[,2]*X2coef.diff[,2]
  EX2.1 <- c(omega1MC.reduced[,1] + omega1MC.reduced[,2:3]%*%x1)
  EX2.2 <- c(omega2MC.reduced[,1] + omega2MC.reduced[,2:3]%*%x1)
  kappa <- X2coef.diff[,1]*EX2.1 + X2coef.diff[,2]*EX2.2 + lambda1
  mopt.s1 <- 0.5*(lambda0 + c(eta1MC.reduced[,4] + eta0MC.reduced[,4])*c(EX2.1) + c(eta1MC.reduced[,5] + eta0MC.reduced[,5])*c(EX2.2) + sqrt(s.x2.sq)*sqrt(2/pi)*exp(-kappa^2/(2*s.x2.sq)) + kappa*(1-2*pnorm(-kappa/sqrt(s.x2.sq))) )
  
  ######Stage 1 when s=-1######
  eta0MC.reduced <- eta0MC[,c(1:3,7:8)]
  omega1MC.reduced <-  omega1MC[,c(1:3)]
  omega2MC.reduced <-  omega2MC[,c(1:3)]
  mopt.sNeg1 <- eta0MC.reduced[,1] + c(eta0MC.reduced[,2:3]%*%x1) + eta0MC.reduced[,4]*(omega1MC.reduced[,1] + c(omega1MC.reduced[,2:3]%*%x1)  ) + eta0MC.reduced[,5]*(omega2MC.reduced[,1] + c(omega2MC.reduced[,2:3]%*%x1)  )
  
  
  optA1 <- getmode(apply(cbind(mopt.sNeg1,mopt.s1),1,which.max))
  optA1 <- c(-1,1)[optA1]
  
  return(list(optA1=optA1))
  
}

ql.glm.ctsX.STARD.CV.v3 = function(dataset,design1Test, design2Test){
  
  n = length(dataset$a1)
  n1 <- sum(dataset$a2); n0 <- n - n1
  z2 <- dataset$z2; z1 <- dataset$z1
  q1 <- ncol(z1); q2 <- ncol(z2); y <- dataset$y
  z2.counterfact <- z2
  z2.counterfact[dataset$a1==1,] <- cbind(z2[dataset$a1==1,9:13],z2[dataset$a1==1,6:8],z2[dataset$a1==1,1:5])
  beta2.hat <- c(solve(t(z2)%*%z2)%*%t(z2)%*%y)
  mean.actual <- c(z2%*%beta2.hat)
  mean.counterfact <- c(z2.counterfact%*%beta2.hat)
  keep.id <- which(mean.actual>=mean.counterfact)
  ymax <- rep(NA,n)
  ymax[keep.id] <- mean.actual[keep.id]
  ymax[-keep.id] <- mean.counterfact[-keep.id]
  beta1.hat <- c(solve(t(z1)%*%z1)%*%t(z1)%*%ymax)
  a1.opt <- 2*as.numeric(c(design1Test%*% beta1.hat[4:6]) >= c(design1Test%*% beta1.hat[1:3])) - 1
  a2.opt <- rep(NA,nrow(design1Test))
  a2.opt[which(a1.opt==-1)] <- -1
  a2.opt[which(a1.opt==1)] <- 2*as.numeric(design2Test[which(a1.opt==1),] %*% beta2.hat[1:5] >= design2Test[which(a1.opt==1),] %*% beta2.hat[9:13]) -1
  
  
  return(list(beta2.hat=beta2.hat, beta1.hat=beta1.hat,a2.opt=a2.opt,a1.opt=a1.opt))
}

bml.glm.ctsX.STARD.CV.v3 = function(dataset,design1Test, design2Test,nsamps=500){
  
  # data
  n = length(dataset$a1)
  n1 <- sum(dataset$a2); n0 <- n - n1
  z2 <- dataset$z2; z1 <- dataset$z1
  q1 <- ncol(z1); q2 <- ncol(z2); y <- dataset$y
  z2.counterfact <- z2
  z2.counterfact[dataset$a1==1,] <- cbind(z2[dataset$a1==1,9:13],z2[dataset$a1==1,6:8],z2[dataset$a1==1,1:5])
  
  X2tX2.inv = solve(t(z2)%*%z2); beta2.hat = X2tX2.inv%*%t(z2)%*%y; X1tX1.inv = solve(t(z1)%*%z1)
  sigma2.sq = 1/rgamma(nsamps,shape=(n-q2)/2,rate=t(y-z2%*%beta2.hat)%*%(y-z2%*%beta2.hat)/2)
  beta2 = t(sapply(1:nsamps, function(samp) mvrnorm(1,beta2.hat,sigma2.sq[samp]*X2tX2.inv)))
  
  sigma1.sq = rep(NA,nsamps); beta1 = matrix(NA,nsamps,q1); 
  ytilde.draw <- rep(NA,n)
  
  for(samp in 1:nsamps){
    
    MEANmat.temp = cbind(c(z2%*%beta2[samp,]), c(z2.counterfact%*%beta2[samp,]))
    stay.id <- as.numeric(MEANmat.temp[,1] >= MEANmat.temp[,2])
    mu.opt <- apply(MEANmat.temp,1,max)
    y2.opt <- rep(0,n)
    y2.opt[stay.id==1] <- y[stay.id==1]
    y2.opt[stay.id==0] <- mu.opt[stay.id==0] + rnorm( n = sum(stay.id==0), sd=sqrt(sigma2.sq[samp]))
    beta1.hat = X1tX1.inv%*%t(z1)%*%y2.opt
    sigma1.sq[samp] = 1/rgamma(1,shape=(n-q1)/2,rate=t(y2.opt-z1%*%beta1.hat)%*%(y2.opt-z1%*%beta1.hat)/2)
    beta1[samp,] = mvrnorm(1,beta1.hat,sigma1.sq[samp]*X1tX1.inv)
    
  }
  
  a1.opt <- apply(cbind(rowMeans(design1Test %*% t(beta1[,(q1/2+1):q1])), rowMeans(design1Test %*% t(beta1[,1:(q1/2)])) ),1,which.max)
  a1.opt <- c(1,-1)[a1.opt]
  a2.opt <- rep(NA,nrow(design2Test))
  a2.opt[a1.opt==-1] <- -1
  a2.opt[a1.opt==1] <- 2*as.numeric(rowMeans(design2Test[a1.opt==1,] %*% t(beta2[,1:5])) >= rowMeans(design2Test[a1.opt==1,] %*% t(beta2[,9:13]))) - 1
  
  return(list( a1.opt=a1.opt, a2.opt=a2.opt, beta2=beta2, beta1=beta1))
}