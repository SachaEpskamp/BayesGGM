#' Bayesian Graphical LASSO.
#'
#' \code{BayesGLASSO} Provides a regularized precision matrix estimate using
#' a Bayesian GLASSO.
#'
#' @param data A data matrix with rows representinh participants and columns
#' representing variables/nodes.
#' @param lambda_shape The rate parameter for the hyperprior of lambda, the parameter
#' that determines the amount of shirnkage.
#' @param lambda_rate The scale parameter for the hyperprior of lambda, the parameter
#' that determines the amount of shirnkage.
#' @param nBurning The number of burn-in iterations for the Gibbs-sampler.
#' @param nIter The number of iterations for the Gibbs-sampler.
#' @param verbose If 'TRUE'displays a progress bar.
#' @param lambda The shrinkage parameter. This value is estimated by default.
#' If given here lambda won't be update as part of the Gibbs-sampler.
#' @return The function returns the following output:
#' \item{Omega}{A dataframe containing the estimated precision matrix for each itteration of the Gibbs-sampler.}
#' \item{pcor}{A dataframe containing the estimated partial correlation matrix for each itteration of the Gibbs-sampler.}
#' \item{lambda_sq}{A dataframe containing the estimates for the local hyper parameter lambda for each itteration of the Gibbs-sampler.}
#' \item{tau_sq}{A dataframe containing the estimates for the global hyper parameter tau for each itteration of the Gibbs-sampler.}
#' \item{Strength}{A dataframe containing the estimated strength centrality for each node for each itteration of the Gibbs-sampler.}
#' \item{Closeness}{A dataframe containing the estimated closeness centrality for each node for each itteration of the Gibbs-sampler.}
#' \item{Betweenness}{A dataframe containing the estimated betweenness centrality for each node for each itteration of the Gibbs-sampler.}
#' \item{optwi}{The point estimate for the precision matrix obtained by taking the mode of the posterior distribution.}
#' \item{optpcor}{The point estimate for the partial correlation matrix obtained by taking the mode of the posterior distribution.}
#' \item{CredInt}{The 95\% Credibility interval for the elements of the precision matrix.}
#' \item{optStrength}{The point estimate for strength centrality of each node based on the point estimate of the partial correlation matrix.}
#' \item{CIStrength}{The 95\% Credibility Interval for the strength centrality of each node.}
#' \item{optCloseness}{The point estimate for closeness centrality of each node based on the point estimate of the partial correlation matrix.}
#' \item{CICloseness}{The 95\% Credibility Interval for the closeness centrality of each node.}
#' \item{optBetween}{The point estimate for betweenness centrality of each node based on the point estimate of the partial correlation matrix.}
#' \item{CIBetween}{The 95\% Credibility Interval for the betweenness centrality of each node.}
#'
#' The following output is only returned in case of missing data.
#'
#' \item{Missing}{Indicator for which observations where missing.}
#' \item{ImputedValues}{The imputed values for each missing datapoint for each itteration of the Gibbs-sampler.}
#' \item{CompleteData}{The original dataset with missing values replaced by the mean imputed value for each missing data point.}
#' \item{OriginalData}{The original dataset with missing values.}
#' @export
BayesGLASSO <- function(
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
    data[,i] <- data[,i] - mean(data[,i],na.rm=T)
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
  Y <- Yorig <- data

  # Indicate missing values in Y
  ind <- which(is.na(Y),arr.ind = T)
  misslabel <- matrix(NA,nrow=nrow(ind),ncol=1)


  # If there are missing values indicate which ones, and intially replace them with the mean score on the variable
  if(nrow(ind)!=0){
    for(i in 1:nrow(ind)){misslabel[i,] <- toString(c(ind[i,1],ind[i,2]))}
  Y[ind] <- colMeans(Y, na.rm=T)[c(as.numeric(ind[,2]))]
  }

  # S:
  S <- t(Y) %*% Y

  # Results:
  Results <- list(
    Omega = array(dim = c(nVar, nVar, nIter)),
    lambda = numeric(nIter),
    pcor =array(dim = c(nVar, nVar, nIter)),
    Strength =matrix(,nrow=nIter,ncol=nVar),
    Closeness =matrix(,nrow=nIter,ncol=nVar),
    Betweenness =matrix(,nrow=nIter,ncol=nVar),
    ImputatedValues = matrix(,nrow=nIter,ncol=nrow(ind))
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

    # Update Y and S

    if(nrow(ind)!=0){
    Ysamp <- mvtnorm::rmvnorm(nObs, sigma=solve(Omega))
    Yitt <- Y
    Yitt[ind] <- Ysamp[ind]
    S <- t(Yitt) %*% Yitt
    Yimp <- matrix(Ysamp[ind],nrow=1)
    }

    # Store:
    if (it > 0){
      Results$Omega[,,it] <- Omega
      Results$lambda[it] <- lambda
      Results$pcor[,,it] <- as.matrix(qgraph:::wi2net(Omega))
      Results$Strength[it,] <- qgraph::centrality(as.matrix(qgraph:::wi2net(Omega)))$OutDegree
      Results$Closeness[it,] <- qgraph::centrality(as.matrix(qgraph:::wi2net(Omega)))$Closeness
      Results$Betweenness[it,] <- qgraph::centrality(as.matrix(qgraph:::wi2net(Omega)))$Betweenness
      if(nrow(ind)!=0){
        Results$ImputatedValues[it,] <- Yimp
      }
    }

    if (verbose){
      setTxtProgressBar(pb, it)
    }
  }

  Results$optwi <- apply(Results$Omega,1:2,function(x)modeest::mlv(x)$M)
  netGLASSO <- qgraph:::wi2net(Results$optwi)

  # Threshold by checking if 0 is in interval:
  SigGLASSO <- apply(Results$Omega,1:2,function(x){
    Quants <- quantile(x,probs = c(0.025,0.975))
    (Quants[1] > 0) | (Quants[2] < 0)
  })

  Results$optpcor <- netGLASSO * SigGLASSO
  Results$CredInt <- apply(Results$Omega,1:2,function(x){
    Quants <- quantile(x,probs = c(0.025,0.975))})

  Results$optStrength <- qgraph::centrality(as.matrix(Results$optpcor))$OutDegree

  # Correct strength estimates per itteration for lack of 0 decision rule by looking at difference in mode of
  # estimates per itteration and the optimal estimates
  Stdiff <- apply(Results$Strength,2,function(x)modeest::mlv(x)$M) - Results$optStrength
  StCorr <- matrix(NA,nrow=nrow(Results$Strength),ncol=ncol(Results$Strength))
  for(i in 1:nrow(Results$Strength)){
    StCorr[i,] <- Results$Strength[i,]- Stdiff
  }
  Results$CIStrength <- apply(StCorr,2,function(x){
    Quants <- quantile(x,probs = c(0.025,0.975))})

  Results$optCloseness <- qgraph::centrality(as.matrix(Results$optpcor))$Closeness

  Cldiff <- apply(Results$Closeness,2,function(x)modeest::mlv(x)$M) - Results$optCloseness
  ClCorr <- matrix(NA,nrow=nrow(Results$Closeness),ncol=ncol(Results$Closeness))
  for(i in 1:nrow(Results$Closeness)){
    ClCorr[i,] <- Results$Closeness[i,]- Cldiff
  }
  Results$CICloseness <- apply(ClCorr,2,function(x){
    Quants <- quantile(x,probs = c(0.025,0.975))})

  Results$optBetweenness <- qgraph::centrality(as.matrix(Results$optpcor))$Betweenness

  Bediff <- apply(Results$Betweenness,2,function(x)modeest::mlv(x)$M) - Results$optBetweenness
  BeCorr <- matrix(NA,nrow=nrow(Results$Betweenness),ncol=ncol(Results$Betweenness))
  for(i in 1:nrow(Results$Betweenness)){
    BeCorr[i,] <- Results$Betweenness[i,]- Bediff
  }
  Results$CIBetweenness <- apply(BeCorr,2,function(x){
    Quants <- quantile(x,probs = c(0.025,0.975))})

  if(nrow(ind)!=0){
  Results$Missings <- ind
  colnames(Results$ImputatedValues) <- misslabel
  Ycomp <- Y
  Ycomp[ind] <- colMeans(Results$ImputatedValues)
  Results$CompleteData <- Ycomp
  Results$OriginalData <- Yorig
  }

  if (verbose) close(pb)

  if(nrow(ind)==0){Results$ImputatedValues <- NULL}

  return(Results)
}
