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
