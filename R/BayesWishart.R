#' Wishart Estimator.
#'
#' \code{BayesWishart} Provides non-regularized precision matrix estimate.
#'
#' @param data A data matrix with rows representinh participants and columns
#' representing variables/nodes.
#' @return Returns dataframes for the precision matrix, partial
#' correlations, and the variability is the estimates of the precision matrix entries
#' with estimates from each iteration of the Gibbs-sampler.
#' @export
BayesWishart <- function(
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
    pcor = as.matrix(qgraph:::wi2net(Omega)),
    Strength <- qgraph::centrality(as.matrix(qgraph:::wi2net(Omega)))$OutDegree,
    Closeness <- qgraph::centrality(as.matrix(qgraph:::wi2net(Omega)))$Closeness,
    Betweenness <- qgraph::centrality(as.matrix(qgraph:::wi2net(Omega)))$Betweenness
  )

  CredInt <- array(NA,dim=c(2,nVar,nVar))
  CredInt[1,,] <- Results$Omega - 1.96*Results$varOmega
  CredInt[2,,] <- Results$Omega + 1.96*Results$varOmega
  Results$CredInt <- CredInt

   return(Results)
}
