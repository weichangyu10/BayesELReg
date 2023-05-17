library(mvtnorm)
library(matrixStats)
library(cvTools)

expit <- function(t){
  
  1/(1+exp(-t))
  
}

softmax <- function(t){
  
  exp(t)/sum(exp(t))
  
  
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

aLC = function (qtau, v, w=NULL){
  if (all(is.na(v))){
    return (list("v"=v, "w"=w));
  }
  if (!is.null(w)){
    ind = setdiff(1:length(v),
                  intersect (which(!is.na(v)), which(!is.na(w))));
    v[ind] = NA;
    w[ind] = NA;
  }
  
  if (all(is.na(v))){
    return (list("v"=v, "w"=w));
  }
  vMax = max (v, na.rm=T);
  v = v - vMax + qtau;
  if (!is.null (w)){
    w[which(v < 0)] = NA;
  }
  v[which(v < 0)] = NA;
  
  if (all (is.na (v))){
    if (is.null(w)){
      return (list ("v"=v, "w"=w));
    }
    return (list("v"=v, "w"=rep(NA, length(v))));
  }
  
  # right align
  nv = length (v);
  vTmp = rep (NA, nv);
  vTmp[(nv + 1 - sum(!is.na(v))):nv] = v[which(!is.na(v))];
  v = vTmp;
  if (!is.null (w)){
    indV = which (!is.na(v));
    indW = which (!is.na(w));
    if (length(indV) != length(indW)){
      print (v);
      print (w);
    }
    wTmp = rep (NA, nv);
    wTmp[indV] = w[indW];
    w = wTmp;
  }
  return (list("v"=v, "w"=w));
}

#Apologies! File "stard_s1.csv" is only available to authorised researchers (with data access rights approved by NIMH-NDA).
#Please email corresponding author with your NIMH-NDA data access certificate to obtain the file.
#Alternatively, you may download the relevant variables directly from NIMH-NDA website.
stard = read.csv ("stard_s1.csv", header=T);
stard = subset(stard, (!is.na(SSRI_l3_ivrc))&(!is.na(SSRI_l2_ivrc)))
n = dim (stard)[1];
#yNms = names(stard)[grep("QSTOT_l2_qs_[0-9]+$", names(stard))]
dts = names(stard)[grep("DATE_l1_qc_[0-9]+$", names(stard))]


l1Q = names(stard)[grep("QSTOT_l1_qs_[0-9]+$", names(stard))];
cgiNms1 = names(stard)[grep("CGI_I_l1_cc_[0-9]+$", names(stard))];
stard$l1MnQ = rep (NA, n); # mean QIDS
stard$l1LstQ = rep (NA, n); # last QIDS
stard$l1FirstQ = rep (NA, n); # First QIDS
stard$dayMax = rep (NA, n);

#W = matrix (NA, nrow=n, ncol=length(cutPts));
for (i in 1:n){
  rawQIDS = as.numeric (stard[i, l1Q]);
  if (any (!is.na(rawQIDS))){
    stard$l1MnQ[i] = mean (rawQIDS, na.rm=T);
    stard$l1LstQ[i] = rawQIDS[max(which(!is.na(rawQIDS)))];
    stard$l1FirstQ[i] = rawQIDS[min(which(!is.na(rawQIDS)))];
    stard$l1dayMax[i] = max (stard[i, dts], na.rm=T);
  }
  rawCGI = as.numeric (stard[i, cgiNms1]);
  if (any (!is.na(rawCGI))){
    stard$l1MnCGI[i] = mean (rawCGI, na.rm=T);
  }
  
}

stard$l1slope <- (stard$l1LstQ - stard$l1FirstQ)/stard$l1dayMax


l2Q = names(stard)[grep("QSTOT_l2_qs_[0-9]+$", names(stard))];
cgiNms2 = names(stard)[grep("CGI_I_l2_cc_[0-9]+$", names(stard))];
dts2 = names(stard)[grep("DATE_l2_qc_[0-9]+$", names(stard))]
stard$l2MnQ = rep (NA, n); # mean QIDS
stard$l2LstQ = rep (NA, n); # last QIDS
stard$l2FirstQ = rep (NA, n); # First QIDS
stard$l2dayMax = rep (NA, n);

#W = matrix (NA, nrow=n, ncol=length(cutPts));
for (i in 1:n){
  rawQIDS = as.numeric (stard[i, l2Q]);
  if (any (!is.na(rawQIDS))){
    stard$l2MnQ[i] = mean (rawQIDS, na.rm=T);
    stard$l2LstQ[i] = rawQIDS[max(which(!is.na(rawQIDS)))];
    stard$l2FirstQ[i] = rawQIDS[min(which(!is.na(rawQIDS)))];
    stard$l2dayMax[i] = max (stard[i, dts2], na.rm=T);
  }
  rawCGI = as.numeric (stard[i, cgiNms2]);
  if (any (!is.na(rawCGI))){
    stard$l2MnCGI[i] = mean (rawCGI, na.rm=T);
  }
  
}

stard$l2slope <- (stard$l2LstQ - stard$l2FirstQ)/stard$l2dayMax


l3Q = names(stard)[grep("QSTOT_l3_qs_[0-9]+$", names(stard))];
stard$l3MnQ = rep (NA, n) # mean QIDS
for (i in 1:n){
  rawQIDS = as.numeric (stard[i, l3Q]);
  if (any (!is.na(rawQIDS))){
    stard$l3MnQ[i] = mean (rawQIDS, na.rm=T);
  }
  
}

A1 = stard$SSRI_l2_ivrc
A2 = stard$SSRI_l3_ivrc
Y = 27 - stard$l3MnQ
Design2 <- cbind(stard$l1MnQ, stard$l1slope,  stard$l2MnQ, stard$l2slope)
Design1 <- cbind(stard$l1MnQ, stard$l1slope)
na.id <- NULL
for(i in 1:n){
  
  if(any(is.na(c(Design2[i,],Y[i])))){
    
    na.id <- c(na.id, i)
    
  }
  
}
na.id <- c(na.id, which((A1==-1)&(A2==1)))

A1 = A1[-na.id]; A2 = A2[-na.id]
Design1 <- cbind(1,Design1[-na.id,])
Design2 <- cbind(1,Design2[-na.id,])
Y <- Y[-na.id]
stage2.MAT <- cbind(apply(Design1,2,function(s){ s*(A2==1)*(A1==1) }), Design2[,4]*(A2==1), Design2[,5]*(A2==1), apply(Design1,2,function(s){ s*(A1==-1)*(A2==-1) }), apply(Design1,2,function(s){ s*(A1==1)*(A2==-1) }), Design2[,4]*(A2==-1), Design2[,5]*(A2==-1) )
stage1.MAT <- cbind(apply(Design1,2,function(s){ s*(A1==-1)}),  apply(Design1,2,function(s){ s*(A1==1)}) )
stage1.MATX2 <- cbind(stage1.MAT,Design2[,4:5])


N <- length(Y)
r1 <- sum((A1==1)&(A2==1))/N
r2 <- sum((A1==1)&(A2==-1))/N
r3 <- sum((A1==-1)&(A2==-1))/N
EV.Proposed <- rep(0,100)
EV.Proposed.Emplik <- rep(0,100)
EV.QL <- rep(0,100)
EV.BML <- rep(0,100)
EV.16 <- rep(0,100)
EV.17 <- rep(0,100)
EV.Neg16 <- rep(0,100)
EV.Neg17 <- rep(0,100)

V=5
for(seed.num in 1:100){
  
  set.seed(seed.num)
  a.a.a1 <- proc.time()[3]
  
  cvSets <- cvFolds(N, V)
  Proposed.num <- rep(0,V); QL.num <- rep(0,V); BML.num <- rep(0,V); Proposed.emplik.num <- rep(0,V)
  Proposed.denom <- rep(0,V); QL.denom <- rep(0,V); BML.denom <- rep(0,V); Proposed.emplik.denom <- rep(0,V)
  
  for(cv in 1:V){
    
    
    testInds <- cvSets$subsets[which(cvSets$which==cv)]
    trainInds <- (1:N)[-testInds]
    
    nTest <- length(testInds)
    nTrain <- length(trainInds)
    
    stage1.MAT.train <- stage1.MAT[trainInds,]
    stage1.MAT.test <- stage1.MAT[testInds,]
    stage2.MAT.train <- stage2.MAT[trainInds,]
    stage2.MAT.test <- stage2.MAT[testInds,]
    Design1.train <- Design1[trainInds,]
    Design1.test <- Design1[testInds,]
    Design2.train <- Design2[trainInds,]
    Design2.test <- Design2[testInds,]
    stage1.MATX2.train <- stage1.MATX2[trainInds,]
    stage1.MATX2.test <- stage1.MATX2[testInds,]
    Y.train <- Y[trainInds]
    Y.test <- Y[testInds]
    A1.train <- A1[trainInds]
    A1.test <- A1[testInds]
    A2.train <- A2[trainInds]
    A2.test <- A2[testInds]
    
    #Fit proposed model
    stage2LM <- lm(Y.train~stage2.MAT.train-1)
    intermediateLM1 <- lm(Design2.train[,4]~stage1.MAT.train-1)
    intermediateLM2 <- lm(Design2.train[,5]~stage1.MAT.train-1)
    eta.coef <- coef(stage2LM)
    omega1.coef <- coef(intermediateLM1)
    omega2.coef <- coef(intermediateLM2)
    etaDraws <- DistributedComp.ctsX.STARD(Y = Y.train, Z2 = stage2.MAT.train, Z1 = stage1.MAT.train, X2 = Design2.train[,4], etaHat = eta.coef, omegaHat = omega1.coef, RSS.Y = sum(summary(stage2LM)$resid^2), RSS.X2 = sum(summary(intermediateLM1)$resid^2), LMdim = ncol(stage2.MAT.train), YorX2 = 1, nsamps = 50000)
    omega1Draws <- DistributedComp.ctsX.STARD(Y = Y.train, Z2 = stage2.MAT.train, Z1 = stage1.MAT.train, X2 = Design2.train[,4], etaHat = eta.coef, omegaHat = omega1.coef, RSS.Y = sum(summary(stage2LM)$resid^2), RSS.X2 = sum(summary(intermediateLM1)$resid^2), LMdim = ncol(stage2.MAT.train), YorX2 = 2, nsamps = 50000)
    omega2Draws <- DistributedComp.ctsX.STARD(Y = Y.train, Z2 = stage2.MAT.train, Z1 = stage1.MAT.train, X2 = Design2.train[,5], etaHat = eta.coef, omegaHat = omega2.coef, RSS.Y = sum(summary(stage2LM)$resid^2), RSS.X2 = sum(summary(intermediateLM2)$resid^2), LMdim = ncol(stage2.MAT.train), YorX2 = 2, nsamps = 50000)
    
    a1.opt.proposed <- rep(NA,nTest)
    for(i in 1:nTest){
      
      a1.opt.proposed[i] <- OptTrt.ctsX.STARDv3(eta1MC = etaDraws$Tsamples[,1:5], eta0MC = etaDraws$Tsamples[,6:13], omega1MC = omega1Draws$Tsamples.X2, omega2MC = omega2Draws$Tsamples.X2, sigma2g1MC = omega1Draws$sigma2g.samples, sigma2g2MC = omega2Draws$sigma2g.samples, x1 = Design1.test[i,-1])$optA1
      
    }
    a2.opt.proposed <- rep(NA,nTest)
    a2.opt.proposed[which(a1.opt.proposed==-1)] <- -1
    a2.opt.proposed[which(a1.opt.proposed==1)] <- 2*as.numeric(apply( (Design2.test[which(a1.opt.proposed==1),]%*%t(etaDraws$Tsamples[,1:5]) >= Design2.test[which(a1.opt.proposed==1),]%*%t(etaDraws$Tsamples[,9:13])), 1, getmode)) - 1
    Proposed.ID <- which((A1.test==a1.opt.proposed)&(A2.test==a2.opt.proposed))
    p.proposed <- rep(NA, length(Proposed.ID))
    if( length((A1.test[Proposed.ID]==1)&(A2.test[Proposed.ID]==1)) > 0 ){
      
      p.proposed[which((A1.test[Proposed.ID]==1)&(A2.test[Proposed.ID]==1))] <- r1
      
    }
    if( length((A1.test[Proposed.ID]==1)&(A2.test[Proposed.ID]==-1)) > 0  ){
      
      p.proposed[which((A1.test[Proposed.ID]==1)&(A2.test[Proposed.ID]==-1))] <- r2
      
    }
    if( length((A1.test[Proposed.ID]==-1)&(A2.test[Proposed.ID]==-1)) > 0 ){
      
      p.proposed[which((A1.test[Proposed.ID]==-1)&(A2.test[Proposed.ID]==-1))] <- r3
      
    }
    Proposed.num[cv] <- sum(Y.test[Proposed.ID]/p.proposed)
    Proposed.denom[cv] <- sum(1/p.proposed)
    
    #Emplik
    start.HMC.time <- Sys.time()
    etaDraws.emplik <- STARD.DrawCoefs.Emplik.HMC(y = Y.train, stagedesign = stage2.MAT.train, estHat = eta.coef, nsamp = 10000)
    omegaDraws1.emplik <- STARD.DrawCoefs.Emplik.HMC(y = Design2.train[,4], stagedesign = stage1.MAT.train, estHat = omega1.coef, nsamp = 10000)
    omegaDraws2.emplik <- STARD.DrawCoefs.Emplik.HMC(y = Design2.train[,5], stagedesign = stage1.MAT.train, estHat = omega2.coef, nsamp = 10000)
    end.HMC.time <- Sys.time()
    
    start.VB.time <- Sys.time()
    etaDraws.emplik <- STARD.DrawCoefs.Emplik.VB(y = Y.train, stagedesign = stage2.MAT.train, estHat = eta.coef, nsamp = 10000)
    omegaDraws1.emplik <- STARD.DrawCoefs.Emplik.VB(y = Design2.train[,4], stagedesign = stage1.MAT.train, estHat = omega1.coef, nsamp = 10000)
    omegaDraws2.emplik <- STARD.DrawCoefs.Emplik.VB(y = Design2.train[,5], stagedesign = stage1.MAT.train, estHat = omega2.coef, nsamp = 10000)
    end.VB.time <- Sys.time()
    
    #omegaDraws.emplik <- STARD.DrawCoefs.X2.Emplik.VB(x21 = Design2.train[,4], x22 = Design2.train[,5], stage1design = stage1.MAT.train, omegaHat1 = omega1.coef, omegaHat2 = omega2.coef, nsamp = 10000)
    ZetaArray <- array(0,dim=c(nrow(omegaDraws.emplik$Tsamples),nTrain,2))
    for(j in 1:nrow(omegaDraws.emplik$Tsamples)){

      ZetaArray[j,,1] <- Design2.train[,4] - c(stage1.MAT.train%*%omegaDraws.emplik$Tsamples[j,1:6])
      ZetaArray[j,,2] <- Design2.train[,5] - c(stage1.MAT.train%*%omegaDraws.emplik$Tsamples[j,7:12])
    }
    a1.opt.proposed.emplik <- rep(NA,nTest)
    for(i in 1:nTest){
      
      a1.opt.proposed.emplik[i] <- OptTrt.STARD.emplik(eta1MC = etaDraws.emplik$Tsamples[,1:5], eta0MC = etaDraws.emplik$Tsamples[,6:13], omega1MC = omegaDraws.emplik$Tsamples[,1:6], omega2MC = omegaDraws.emplik$Tsamples[,7:12], x1 = Design1.test[1,2:3], XiMatrix1 = ZetaArray[,,1], XiMatrix2 = ZetaArray[,,2])$optA1
      
    }
    a2.opt.proposed.emplik <- rep(NA,nTest)
    a2.opt.proposed.emplik[which(a1.opt.proposed.emplik==-1)] <- -1
    a2.opt.proposed.emplik[which(a1.opt.proposed.emplik==1)] <- 2*as.numeric(apply( (Design2.test[which(a1.opt.proposed.emplik==1),]%*%t(etaDraws.emplik$Tsamples[,1:5]) >= Design2.test[which(a1.opt.proposed.emplik==1),]%*%t(etaDraws.emplik$Tsamples[,9:13])), 1, getmode)) - 1
    Proposed.Emplik.ID <- which((A1.test==a1.opt.proposed.emplik)&(A2.test==a2.opt.proposed.emplik))
    p.proposed.Emplik <- rep(NA, length(Proposed.Emplik.ID))
    if( length((A1.test[Proposed.Emplik.ID]==1)&(A2.test[Proposed.Emplik.ID]==1)) > 0 ){
      
      p.proposed.Emplik[which((A1.test[Proposed.Emplik.ID]==1)&(A2.test[Proposed.Emplik.ID]==1))] <- r1
      
    }
    if( length((A1.test[Proposed.Emplik.ID]==1)&(A2.test[Proposed.Emplik.ID]==-1)) > 0  ){
      
      p.proposed.Emplik[which((A1.test[Proposed.Emplik.ID]==1)&(A2.test[Proposed.Emplik.ID]==-1))] <- r2
      
    }
    if( length((A1.test[Proposed.Emplik.ID]==-1)&(A2.test[Proposed.Emplik.ID]==-1)) > 0 ){
      
      p.proposed.Emplik[which((A1.test[Proposed.Emplik.ID]==-1)&(A2.test[Proposed.Emplik.ID]==-1))] <- r3
      
    }
    Proposed.emplik.num[cv] <- sum(Y.test[Proposed.Emplik.ID]/p.proposed.Emplik)
    Proposed.emplik.denom[cv] <- sum(1/p.proposed.Emplik)

    #Q-Learning
    QLobj <- ql.glm.ctsX.STARD.CV.v3(dataset = list(a1=A1.train, a2 = A2.train, design1 = Design1.train, z2 = stage2.MAT.train, z1=stage1.MAT.train, y = Y.train), design1Test = Design1.test, design2Test = Design2.test)
    #QLobj <- ql.glm.ctsX.STARD.CV.v6(dataset = list(a1=A1.train, a2 = A2.train, design1 = Design1.train, z2 = stage2.MAT.train, z1=stage1.MAT.train, y = Y.train), design1Test = Design1.test, design2A1Test = stage1.MATX2.test)
    QL.ID <- which((A1.test==QLobj$a1.opt)&(A2.test==QLobj$a2.opt))
    p.QL <- rep(NA, length(QL.ID))
    if( length((A1.test[QL.ID]==1)&(A2.test[QL.ID]==1)) > 0 ){
      
      p.QL[which((A1.test[QL.ID]==1)&(A2.test[QL.ID]==1))] <- r1
      
    }
    if( length((A1.test[QL.ID]==1)&(A2.test[QL.ID]==-1)) > 0  ){
      
      p.QL[which((A1.test[QL.ID]==1)&(A2.test[QL.ID]==-1))] <- r2
      
    }
    if( length((A1.test[QL.ID]==-1)&(A2.test[QL.ID]==-1)) > 0 ){
      
      p.QL[which((A1.test[QL.ID]==-1)&(A2.test[QL.ID]==-1))] <- r3
      
    }
    QL.num[cv] <- sum(Y.test[QL.ID]/(p.QL))
    QL.denom[cv] <- sum(1/(p.QL))
    
    
    #BML
    BMLobj <- bml.glm.ctsX.STARD.CV.v3(dataset = list(a1=A1.train, a2 = A2.train, design1 = Design1.train, z2 = stage2.MAT.train, z1=stage1.MAT.train, y = Y.train), design1Test = Design1.test, design2Test = Design2.test, nsamps = 50000)
    BML.ID <- which((A1.test==BMLobj$a1.opt)&(A2.test==BMLobj$a2.opt))
    p.BML <- rep(NA, length(BML.ID))
    if( length((A1.test[BML.ID]==1)&(A2.test[BML.ID]==1)) > 0 ){
      
      p.BML[which((A1.test[BML.ID]==1)&(A2.test[BML.ID]==1))] <- r1
      
    }
    if( length((A1.test[BML.ID]==1)&(A2.test[BML.ID]==-1)) > 0  ){
      
      p.BML[which((A1.test[BML.ID]==1)&(A2.test[BML.ID]==-1))] <- r2
      
    }
    if( length((A1.test[BML.ID]==-1)&(A2.test[BML.ID]==-1)) > 0 ){
      
      p.BML[which((A1.test[BML.ID]==-1)&(A2.test[BML.ID]==-1))] <- r3
      
    }
    BML.num[cv] <- sum(Y.test[BML.ID]/(p.BML))
    BML.denom[cv] <- sum(1/(p.BML))
    
      
      
    }
    
    
  EV.Proposed[seed.num] <- sum(Proposed.num)/sum(Proposed.denom)
  EV.QL[seed.num] <- sum(QL.num)/sum(QL.denom)
  EV.Proposed.Emplik[seed.num] <- sum(Proposed.emplik.num)/sum(Proposed.emplik.denom)
  EV.BML[seed.num] <- sum(BML.num)/sum(BML.denom)
  
  
  if(seed.num%%10==0){
    
    save.image("STARDanalysisWithOutliersv3WithCVandEmplik.Rdata")
    
  }
  
  
}


