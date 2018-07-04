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
