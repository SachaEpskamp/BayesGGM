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
  lambda_shape = 1, # Hyperprior lambda  rate??
  lambda_rate = 0.01, # Hyperprior lambda scale??
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
 # ????
  
  # Tau??
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
     
     gamma <- rgamma(1, shape = nObs / 2 + 1, rate = (S22 + lambda) / 2)
     C <- solve((S22 + lambda) * solve(Omega11)+ solve(diag(c(Tau12))))
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