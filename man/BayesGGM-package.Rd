\name{BayesGGM-package}
\alias{BayesGGM-package}
\alias{BayesGGM}
\docType{package}
\title{
  Bayesian Estimators for Gaussian Graphical Models
}
\description{
  BayesGGM is a R package to estimate Gaussian Graphical Models (GGMs) using 6 different Bayesian estimation methods.
  Researchers can choose between 1) a Bayesian GLASSO, 2) an Adaptive Bayesian GLASSO, 3) a Graphical Horseshoe, 4)
  a (Ridge-type) penalized estimator based on an Eigenvalue decompostion of the Wishart prior, 5) an adaptive (Ridge-type) penalized estimator based on an Eigenvalue decomposition of the Wishart prior, and 6) a normal Wishart prior with no regularization.
  In addition, BayesGGM also provide 95\% Credibility Intervals for the strength, closeness, and betweenness centrality, and imputes missing data.
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Joran Jongerling.

Maintainer: Joran Jongerling <jongerling@essb.eur.nl>
}
\references{
  Main literature:

  Kuismin, M. & Sillanpaa, M. (2016). Use of Wishart Prior and Simple Extensions for Sparse Precision Matrix Estimation.

  Li, Y. (2018). The Graphical Horeshoe Estimator for Inverse Coavariance Matrices.

  Wang, H. (2012). Bayesian graphical lasso models and efficient posterior computation. Bayesian Analysis, 7(4), 867-886.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.
  }
}
