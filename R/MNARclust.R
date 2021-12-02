##' MNARclust.
##'
##'
##'
##' \tabular{ll}{
##'   Package: \tab MNARclust\cr
##'   Type: \tab Package\cr
##'   Version: \tab 1.1.0\cr
##'   Date: \tab 2021-12-01\cr
##'   License: \tab GPL-3\cr
##'   LazyLoad: \tab yes\cr
##' }
##'
##' @description
##' Clustering method to analyze continuous or mixed-type data with missingness. The missingness mechanism can be non ignorable. The approach considers a semi-parametric mixture model.
##'
##' @name MNARclust-package
##' @aliases MNARclust
##' @rdname MNARclust-package
##' @docType package
##' @keywords package
##' @import parallel
##' @import Rcpp
##' @importFrom sn rsn dsn
##' @importFrom rmutil rlaplace dlaplace
##' @importFrom stats dnorm dt integrate na.omit rnorm rt runif sd
##' @useDynLib  MNARclust
NULL

##' Echocardiogram data set
##'
##'
##' @description All the patients suffered heart attacks at some point in the past. Some are still alive and some are not. The survival and still-alive variables, when taken together, indicate whether a patient survived for at least one year following the heart attack.
##' @details This data set arise from the UCI machine learning repository (more details on this data set are presented http://archive.ics.uci.edu/ml/datasets/Echocardiogram)
##' @references  Salzberg, S. (1988). Exemplar-based learning: Theory and implementation (Technical Report TR-10-88). Harvard University, Center for Research in Computing Technology, Aiken Computation Laboratory (33 Oxford Street; Cambridge, MA 02138).
##' @name echo
##' @format A data frame with 132 observations on 13 variables
##' (more details on this data set are presented in http://archive.ics.uci.edu/ml/datasets/Echocardiogram).
##'
##' @docType data
##' @keywords datasets
##'
##' @examples
##' data(echo)
NULL

##' Clustering function
##'
##' @description
##' Clustering method to analyze continuous or mixed-type data with missingness. The missingness mechanism can be non ignorable. The approach considers a semi-parametric mixture model.
##'
##' @param x matrix used for clustering
##' @param K number of components
##' @param nbinit number of random starting points
##' @param nbCPU number of CPU used for parallel computing (only Unix and Linux systems are allowed)
##' @param tol stopping rule
##' @param band bandwidth (numeric vector).
##' @param seedvalue value of the seed (used to set the initializations of the MM algorithm)
##'
##' @return Returns a list containing the proportions (proportions), matrix of probabilities of missngness (rho), the posterior probabilities of classification (classproba), the partition (zhat) and the logarithme of the smoothed-likelihood (logSmoothlike)
##' @references Clustering Data with Non-Ignorable Missingness using Semi-Parametric Mixture Models, Marie Du Roy de Chaumaray and  Matthieu Marbac <arXiv:2009.07662>.
##'
##' @examples
##' \donttest{
##' set.seed(123)
##' # Data generation
##' ech <- rMNAR(n=100, K=2, d=4, delta=2, gamma=2)
##' # Clustering
##' res <- MNARcluster(ech$x, K=2)
##' # Confusion matrix between the estimated and the true partiion
##' table(res$zhat, ech$z)
##' }
##' @export
MNARcluster <- function(x, K, nbinit = 20, nbCPU = 1, tol = 0.01, band = band.default(x), seedvalue=123){
  ech <- builddata(x)
  all.res <- list()
  set.seed(seedvalue)
  if (length(ech$contvbles) == ech$d){
    all.init <- replicate(nbinit, start.paramCont(ech, K, band), simplify = FALSE)
    if (Sys.info()["sysname"] != "Windows"){
      all.res <- mclapply(all.init, algo_NP_onealgoCont, ech = ech, tol = tol, mc.cores = nbCPU)
    }else{
      all.res <- lapply(all.init, algo_NP_onealgoCont, ech = ech, tol = tol)
    }
  }else if (length(ech$catvbles) == ech$d){
    all.init <- replicate(nbinit, start.paramCat(ech, K, band), simplify = FALSE)
    if (Sys.info()["sysname"] != "Windows"){
      all.res <- mclapply(all.init, algo_NP_onealgoCat, ech = ech, tol = tol,  mc.cores = nbCPU)
    }else{
      all.res <- lapply(all.init, algo_NP_onealgoCat, ech = ech, tol = tol)
    }
  }else{
    all.init <- replicate(nbinit, start.paramMixed(ech, K, band), simplify = FALSE)
    if (Sys.info()["sysname"] != "Windows"){
      all.res <- mclapply(all.init, algo_NP_onealgoMixed, ech = ech, tol = tol, mc.cores = nbCPU)
    }else{
      all.res <- lapply(all.init, algo_NP_onealgoMixed, ech = ech, tol = tol)
    }
  }
  correct <- mean(sapply(all.res, function(u) u$prec) <= sapply(all.res, function(u) u$logSmoothlike))
  all.res <- all.res[[which.max(sapply(all.res, function(u) u$logSmoothlike))]]
  all.res$correct <- correct
  out <- list(proportions = all.res$param$pi,
              rho = all.res$param$rho,
              classproba = all.res$param$weights,
              zhat = all.res$zhat,
              logSmoothlike = all.res$logSmoothlike)
}
