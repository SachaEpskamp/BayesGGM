#' Regularized Wishart Estimator.
#'
#' \code{BayesRegularizedWishart} Provides a ridge type regularized precision matrix estimate using
#' an eigenvalue and eigenvector decomposition of the sample covariance matrix.
#'
#' @param data A data matrix with rows representinh participants and columns
#' representing variables/nodes.
#' @param adaptive If 'TRUE' chose the amount of shrinkage based on the ration of n to p.
#' @param adfact the ration between the two shrinkage hyper-parameters used.
#' @return Returns dataframes for the precision matrix, partial
#' correlations, and an overview of the edges fixed to 0 with estimates from each
#' iteration of the Gibbs-sampler.
#' @export
BayesRegularizedWishart <- function(
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
    pcor = as.matrix(qgraph:::wi2net(OmegaMAPSparse)),
    Strength <- qgraph::centrality(as.matrix(qgraph:::wi2net(OmegaMAPSparse)))$OutDegree,
    Closeness <- qgraph::centrality(as.matrix(qgraph:::wi2net(OmegaMAPSparse)))$Closeness,
    Betweenness <- qgraph::centrality(as.matrix(qgraph:::wi2net(OmegaMAPSparse)))$Betweenness
  )


  return(Results)
}








