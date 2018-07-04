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
