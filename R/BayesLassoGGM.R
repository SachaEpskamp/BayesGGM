sampleTau <- function(Omega, lambda){
  nVar <- ncol(Omega)
  Tau <- matrix(,nVar,nVar)
  for (i in 1:nVar){
    u <- statmod::rinvgauss(1,sqrt(lambda^2/Omega[i,i]^2),lambda^2)
    Tau[i,i] <- 1/u
    for (j in 1:nVar){
      u <- statmod::rinvgauss(1,sqrt(lambda^2/Omega[i,j]^2),lambda^2)
      Tau[i,j] <- Tau[j,i] <- 1/u
    }
  }
  Tau
}

BayesLassoGGM <- function(
  data,
  lambda_shape = 1, # Hyperprior lambda  rate
  lambda_rate = 0.01, # Hyperprior lambda scale
  nBurnin = 10000,
  nIter = 10000,
  verbose = TRUE,
  lambda # if give, do not update
){
  # Main literature:
  # Wang, H. (2012). Bayesian graphical lasso models and efficient posterior computation. Bayesian Analysis, 7(4), 867-886.
  
  # Some parameters:
  nVar <- ncol(data)
  nObs <- nrow(data)
  
  # Center data:
  for (i in 1:nVar){
    data[,i] <- data[,i] - mean(data[,i])
  }
  
  # Use glasso with some lambda to get starting values:
  glassoRes <- glasso::glasso(cov(data,use="pairwise.complete.obs"),rho=0.1)
  
  # Starting values:
  # Inverse:
  Omega <- glassoRes$wi
  
  # Lasso parameter:
  sampleLambda <- missing(lambda)
  
  if (sampleLambda){
    lambda <- rgamma(1,shape = lambda_shape,rate = lambda_rate)    
  }
 
  # Tau
  Tau <- sampleTau(Omega,lambda)
  
  # Y:
  Y <- data
  
  # S:
  S <- t(Y) %*% Y
  
  # Results:
  Results <- list(
    Omega = array(dim = c(nVar, nVar, nIter)),
    lambda = numeric(nIter),
    pcor =array(dim = c(nVar, nVar, nIter))
  )
  
  # Progress bar:
  if (verbose){
    pb <- txtProgressBar(min = -nBurnin, max = nIter, style = 3,initial = 0)  
  }
  
  # Start sampling:
  for (it in (-nBurnin):nIter){
    # For every node:
    for (i in 1:nVar){
      # Partition matrices:
      Omega11 <- Omega[-i,-i]
      Omega12 <- Omega[i,-i]
      Omega21 <- Omega[-i,i]
      Omega22 <- Omega[i,i]
      
      S11 <- S[-i,-i]
      S12 <- S[i,-i]
      S21 <- S[-i,i]
      S22 <- S[i,i]
      
      Tau11 <- Tau[-i,-i]
      diag(Tau11) <- 0
      Tau12 <- Tau[i,-i]
      Tau21 <- Tau[-i,i]
      Tau22 <- 0
      
      gamma <- rgamma(1, shape = nObs/2 + 1, rate = (S22 + lambda) / 2)
      C <- solve((S22 + lambda) * solve(Omega11)+ solve(diag(c(Tau12))))
      #C <- solve((S22 + lambda) * solve(Omega11)+ diag(1/c(Tau12)))
      C <- (C + t(C))/2
      beta <- mvtnorm::rmvnorm(1,-C %*% S21, C)
      Omega[-i,i] <- Omega[i,-i] <- beta
      Omega[i,i] <- gamma + t(t(beta)) %*% solve(Omega11) %*% t(beta)
    }
    
    # Update tau:
    Tau <- sampleTau(Omega,lambda)
    
    # Update lambda:
    if (sampleLambda){
      lambda <- rgamma(1, shape =lambda_shape + nVar * (nVar+1) / 2, rate = lambda_rate + sum(abs(Omega))/2)
     }
    
    # Store:
    if (it > 0){
      Results$Omega[,,it] <- Omega
      Results$lambda[it] <- lambda
      Results$pcor[,,it] <- as.matrix(qgraph:::wi2net(Omega))
    }
    
    if (verbose){
      setTxtProgressBar(pb, it)
    }
  }
  if (verbose) close(pb)
  
  return(Results)
}

sampleTauAdaptive <- function(Omega, lambda, Cadj){
  nVar <- ncol(Omega)
  Tau <- matrix(0,nVar,nVar)
  for (i in 1:nVar){
    u <- statmod::rinvgauss(1,min(sqrt(lambda[i,i]^2/Cadj[i,i]^2),10^12),lambda[i,i]^2)
    Tau[i,i] <- max(10^-12,1/u)
    
    for (j in 1:nVar){
      u <- statmod::rinvgauss(1,min(sqrt(lambda[i,j]^2/Cadj[i,j]^2),10^12),lambda[i,j]^2)
      Tau[i,j] <- Tau[j,i] <- max(10^-12,1/u)
    }
  }
  Tau
}

BayesAdaptiveLassoGGM <- function(
  data,
  lambda_shape = 0.01, # Hyperprior lambda  rate
  lambda_rate = 0.000001, # Hyperprior lambda scale
  nBurnin = 10000,
  nIter = 10000,
  verbose = TRUE,
  lambdaDiag=1
){
  # Main literature:
  # Wang, H. (2012). Bayesian graphical lasso models and efficient posterior computation. Bayesian Analysis, 7(4), 867-886.
  
  # Some parameters:
  nVar <- ncol(data)
  nObs <- nrow(data)
  
  # Center data:
  for (i in 1:nVar){
    data[,i] <- data[,i] - mean(data[,i])
  }
  
  # Use glasso with some lambda to get starting values:
  glassoRes <- glasso::glasso(cov(data,use="pairwise.complete.obs"),rho=0.1)
  
  # Starting values:
  # Inverse:
  Omega <- glassoRes$wi
  
  lambda <- matrix(NA,nrow(Omega),nrow(Omega))
  Cadj <- matrix(NA,nrow(Omega),nrow(Omega))
  
  
  
  for (i in 1:nrow(Omega)){
    for (j in 1:nrow(Omega)){
      Cadj[i,j] <-  max(abs(Omega[i,j]),10^-12)
      lambda[i,j] <- lambda[j,i] <- rgamma(1,shape = 1 + lambda_shape,rate = Cadj[i,j] + lambda_rate) 
    }
  }
  diag(lambda) <- lambdaDiag <- 1
  
  # Tau??
  Tau <- sampleTauAdaptive(Omega,lambda, Cadj)
  
  # Y:
  Y <- data
  
  # S:
  S <- t(Y) %*% Y
  
  # Results:
  Results <- list(
    Omega = array(dim = c(nVar, nVar, nIter)),
    lambda = array(dim = c(nVar, nVar, nIter)),
    pcor =array(dim = c(nVar, nVar, nIter))
  )
  
  
  # Progress bar:
  if (verbose){
    pb <- txtProgressBar(min = -nBurnin, max = nIter, style = 3,initial = 0)  
  }
  
  # Start sampling:
  for (it in (-nBurnin):nIter){
    
    # For every node:
    for (i in 1:nVar){
      # Partition matrices:
      Omega11 <- Omega[-i,-i]
      Omega12 <- Omega[i,-i]
      Omega21 <- Omega[-i,i]
      Omega22 <- Omega[i,i]
      
      S11 <- S[-i,-i]
      S12 <- S[i,-i]
      S21 <- S[-i,i]
      S22 <- S[i,i]
      
      Tau11 <- Tau[-i,-i]
      diag(Tau11) <- 0
      Tau12 <- Tau[i,-i]
      Tau21 <- Tau[-i,i]
      Tau22 <- 0
      
      gamma <- rgamma(1, shape = nObs / 2 + 1, rate = (S22 + lambdaDiag) / 2)
      C <- solve((S22 + lambdaDiag) * solve(Omega11) + diag(c(1/Tau12)))
      C <- (C + t(C))/2
      beta <- mvtnorm::rmvnorm(1,-C %*% S21, C)
      Omega[-i,i] <- Omega[i,-i] <- beta
      Omega[i,i] <- gamma + t(t(beta)) %*% solve(Omega11) %*% t(beta)
      
    }
    
    # Update lambda and tau:
    for (i in 1:nrow(Omega)){
      for (j in 1:nrow(Omega)){
        Cadj[i,j] <-  max(abs(Omega[i,j]),10^-12)
        lambda[i,j] <- lambda[j,i] <- rgamma(1,shape = 1 + lambda_shape,rate = Cadj[i,j] + lambda_rate) 
      }
    }
    diag(lambda) <- lambdaDiag <- 1
    
    # Tau??
    Tau <- sampleTauAdaptive(Omega,lambda, Cadj)
    
    # Store:
    if (it > 0){
      Results$Omega[,,it] <- Omega
      Results$lambda[,,it] <- lambda
      Results$pcor[,,it] <- as.matrix(qgraph:::wi2net(Omega))
    }
    
    if (verbose){
      setTxtProgressBar(pb, it)
    }
  }
  if (verbose) close(pb)
  
  return(Results)
}

BayesHoreshoeGGM <- function(
  data,
  nBurnin = 10000,
  nIter = 10000,
  verbose = TRUE
){
  # Main literature:
  # Li, Y. (2018). The Graphical Horeshoe Estimator for Inverse Coavariance Matrices.
  
  # Some parameters:
  nVar <- ncol(data)
  nObs <- nrow(data)
  
  # Center data:
  for (i in 1:nVar){
    data[,i] <- data[,i] - mean(data[,i])
  }
  
  ## Use glasso with some lambda to get starting values:
  glassoRes <- glasso::glasso(cov(data,use="pairwise.complete.obs"),rho=0.1)
  
  # Starting values:
  # Inverse:
  Omega <- glassoRes$wi
  
  # Tuning and augmentation parameters
  lambda_sq <- matrix(1,nVar,nVar)
  nu <- matrix(1,nVar,nVar)
  tau_sq <- 1
  xi <- 1
  
  # Y:
  Y <- data
  
  # S
  S <- t(Y) %*% Y
  
  # Results:
  Results <- list(
    Omega = array(dim = c(nVar, nVar, nIter)),
    lambda_sq = array(dim = c(nVar, nVar, nIter)),
    tau_sq = numeric(nIter),
    pcor =array(dim = c(nVar, nVar, nIter))
  )
  
  # Progress bar:
  if (verbose){
    pb <- txtProgressBar(min = -nBurnin, max = nIter, style = 3,initial = 0)  
  }
  
  # Start sampling:
  for (it in (-nBurnin):nIter){
    # For every node:
    for (i in 1:nVar){
      # Partition matrices:
      Omega11 <- Omega[-i,-i]
      Omega12 <- Omega[i,-i]
      Omega21 <- Omega[-i,i]
      Omega22 <- Omega[i,i]
      
      S11 <- S[-i,-i]
      S12 <- S[i,-i]
      S21 <- S[-i,i]
      S22 <- S[i,i]
      
      lambda_sq_12 <- lambda_sq[-i,i]
      nu_12 <- nu[-i,i]
      
      gamma <- rgamma(1, shape = (nObs/2 + 1), rate = S22/2)
      C <- solve((S22) * solve(Omega11) + diag(1/c(lambda_sq_12*tau_sq)))
      C <- (C + t(C))/2
      beta <- mvtnorm::rmvnorm(1,-C %*% S21, C)
      Omega[-i,i] <- Omega[i,-i] <- beta
      Omega[i,i] <- gamma + t(t(beta)) %*% solve(Omega11) %*% t(beta)
      rate <- Omega[-i,i]^2/(2*tau_sq) + 1/nu_12
      lambda_sq_12 <- 1/rgamma((nVar-1),shape=1, rate=rate)
      nu_12 <- 1/rgamma((nVar-1),shape=1, rate=1+1/lambda_sq_12)
      lambda_sq[-i,i] <- lambda_sq[i,-i] <- lambda_sq_12
      nu[-i,i] <- nu[i,-i] <- nu_12
    }
    
    # Update tau_sq and xi:
    Omega.vec <- Omega[lower.tri(Omega)]
    lambda.sq.vec <- lambda_sq[lower.tri(lambda_sq)]
    rate <- 1/xi + sum(Omega.vec^2/(2*lambda.sq.vec))
    tau_sq <- 1/rgamma(1, shape = (nVar * (nVar-1)/2 + 1)/2, rate=rate)
    xi <- 1/rgamma(1, shape=1, rate=1+1/tau_sq)
    
    # Store:
    if (it > 0){
      Results$Omega[,,it] <- Omega
      Results$lambda_sq[,,it] <- lambda_sq
      Results$pcor[,,it] <- as.matrix(qgraph:::wi2net(Omega))
      Results$tau_sq[it] <- tau_sq
    }
    
    if (verbose){
      setTxtProgressBar(pb, it)
    }
  }
  if (verbose) close(pb)
  
  return(Results)
}

BayesWishGGM <- function(
  data,
  mu0=0
){
  
  # Some parameters:
  nVar <- ncol(data)
  nObs <- nrow(data)
  
  # Center data:
  for (i in 1:nVar){
    data[,i] <- data[,i] - mean(data[,i])
  }
  
  # Determine sample means
  meanVar <- rep(NA, nVar)
  for (i in 1:nVar){
    meanVar <- mean(data[,i])
  }
  
  
  # Use glasso with some lambda to get starting values:
  glassoRes <- glasso::glasso(cov(data,use="pairwise.complete.obs"),rho=0.1)
  
  # Starting values:
  # Inverse:
  Omega <- glassoRes$wi
  
  # Y:
  Y <- data
  
  # S:
  S = crossprod(Y)
  
  # mu
  mu <- rep(mu0, nVar)
  diff <- meanVar-mu
  
  # Prior scale and df for wishart
  Scale <- Omega
  df0 <- nVar
  
  # Posterior scale and df for Wishart
  ScaleUpd <- solve(nObs*S + solve(Scale))
  df <- df0 + nObs
  
  
  Omega <- (df+nObs-nVar-1)*ScaleUpd
  varOmega <- matrix(NA,nrow=nVar, ncol=nVar)
  
  for (i in 1:nVar){
    for (j in 1:nVar){
      varOmega[i,j] <- (df+nObs)*(Omega[i,j])^2 + Omega[i,i] + Omega[j,j]
    }
  }
  
  # Results:
  Results <- list(
    Omega = Omega,
    varOmega = varOmega,
    pcor = as.matrix(qgraph:::wi2net(Omega))
  )
  
  
  return(Results)
}

CleanUpEBIC = function(ThetaPost,Y,gamma=NULL){
  
  if(is.null(gamma)) gamma = .5
  
  p = ncol(ThetaPost)
  n = nrow(Y)
  
  NotPosDef = 0 # An indicator indicating if we have to do an ad hoc correction to the final estimate to 
  #guarantee positive definiteness of the estimate
  
  S = crossprod(Y)/n
  
  nmb = p*(p-1)/2
  
  A = abs(ThetaPost)
  
  Ak = A[upper.tri(A)]
  
  Ak = sort(Ak)
  
  NmbOfNonZeros = rep(0,nmb)
  
  EBIC = rep(0,nmb)
  
  # This is a very greedy procedure and should only be used in moderate size (say, p = 200)
  
  for(i in 1:nmb){
    
    ThetaPost.threshold = ifelse(A <= Ak[i], 0, ThetaPost)
    
    diag(ThetaPost.threshold) = diag(ThetaPost)
    
    # R seems to have hard time to make symmetric matrices symmetric:
    
    if(!isSymmetric(ThetaPost.threshold)) ThetaPost.threshold = .5*(ThetaPost.threshold + t(ThetaPost.threshold))
    
    # In the case of non positive definite MAP:
    
    if(min(eigen(ThetaPost.threshold)$values) <= 0) ThetaPost.threshold = ThetaPost.threshold + (abs(min(eigen(ThetaPost.threshold)$values)) + 0.01)*diag(1,p)
    
    log.like = (n/2)*(log(det(ThetaPost.threshold)) - sum(diag(S%*%ThetaPost.threshold)))
    
    E = ifelse(ThetaPost.threshold != 0, 1, 0); diag(E) = 0
    
    E.card = sum(E[upper.tri(E)])
    
    NmbOfNonZeros[i] = E.card
    
    EBIC[i] = -2*log.like + E.card*log(n) + 4*E.card*gamma*log(p)
    
  }
  
  ThetaEst = matrix(0,p,p)
  
  a = which.min(EBIC)
  
  ThetaEst = ifelse(abs(ThetaPost) > A[a], ThetaPost, 0); diag(ThetaEst) = diag(ThetaPost)
  
  if(!isSymmetric(ThetaEst)) ThetaEst = .5*(ThetaEst + t(ThetaEst))
  
  # A small correction in the case of not positive definite estimate
  
  if(min(eigen(ThetaEst)$values) < 0){
    
    ThetaEst = ThetaEst + (abs(min(eigen(ThetaEst)$values)) + 0.01)*diag(1,p)
    
    NotPosDef = 1
    
  }
  
  results = list("EBIC" = EBIC, "Theta" = ThetaEst, "Zeros" = NmbOfNonZeros, "NotPosDef" = NotPosDef)
  
  return(results)
  
}

BayesWishEVGGM <- function(
  data,
  adaptive=FALSE,
  adfact=2
){
  
  # Some parameters:
  nVar <- ncol(data)
  nObs <- nrow(data)
  
  # Center data:
  for (i in 1:nVar){
    data[,i] <- data[,i] - mean(data[,i])
  }
  
  
  # Use glasso with some lambda to get starting values:
  glassoRes <- glasso::glasso(cov(data,use="pairwise.complete.obs"),rho=0.1)
  
  # Starting values:
  # Inverse:
  Omega <- glassoRes$wi
  
  # Y:
  Y <- data
  
  S = crossprod(Y)/nObs
  
  ConditionNumberMLE = max(eigen(S)$values)/min(eigen(S)$values)
  
  I = diag(1,nVar)
  
  r = nVar
  
  alpha = r
  
  if(adaptive==TRUE){
    
    if(nObs<nVar){alpha <- adfact*r}
    else{alpha <-r-nVar-1}
  }
  
  
  OmegaMAP = (r + nObs - nVar - 1)*solve(alpha*I + nObs*S)
  
  OmegaMAP = as.matrix(OmegaMAP)
  
  condnmb = max(eigen(OmegaMAP)$values)/min(eigen(OmegaMAP)$values)
  
  b = condnmb
  
  r = r + 0.05
  alpha = r
  
  if(adaptive==TRUE){
    
    if(nObs<nVar){alpha <- adfact*r}
    else{alpha <-r-nVar-1}
  }
  
  while(condnmb > b - 2/sqrt(nObs)){
    
    OmegaMAP = (r + nObs - nVar - 1)*solve(alpha*I + nObs*S)
    
    condnmb = max(eigen(OmegaMAP)$values)/min(eigen(OmegaMAP)$values)
    
    r = r + 0.5
    alpha = r
    
    if(adaptive==TRUE){
      
      if(nObs<nVar){alpha <- adfact*r}
      else{alpha <-r-nVar-1}
    }
    
  }
  
  NewConditionNumberMAP = max(eigen(OmegaMAP)$values)/min(eigen(OmegaMAP)$values)
  
  # Produce the sparse network. We use the same value for the EBIC in our procedure and in both glasso and
  # Meinshausen & Buhlmann approximation
  
  gamma = 0.6
  
  Results = CleanUpEBIC(OmegaMAP,Y,gamma=gamma)
  
  AdHocFixedEstimates = Results$NotPosDef
  
  OmegaMAPSparse = Results$Theta
  
  ConditionNumberMAPSparse = max(eigen(OmegaMAPSparse)$values)/min(eigen(OmegaMAPSparse)$values)
  
  
  # Results:
  Results <- list(
    Omega = OmegaMAPSparse,
    AdHocFixedEstimates = AdHocFixedEstimates,
    pcor = as.matrix(qgraph:::wi2net(OmegaMAPSparse))
  )
  
  
  return(Results)
}