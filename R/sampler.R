##' Function used to simulate data from mixture model with specific missingness mechanism
##'
##' @description
##' Generation of data set to perform the simulation presented in Section 4.1 of Du Roy de Chaumaray (2020)
##'
##' @param n sample size (numeric of length 1)
##' @param K number of clusters (numeric of length 1)
##' @param d number of variables (numeric of length 1)
##' @param delta tuning parameter to define the rate of misclassification (numeric of length 1)
##' @param gamma tuning parameter to define the rate of missingness (numeric of length 1)
##' @param law specifies the distribution of the variables within components (character that must be equal to gauss, student, laplace or skewgauss)
##' @param linkmissing specify the missingness mechanism (character that must be equal to MCAR, logit-Z, logit-X or censoring)
##'
##' @return rMNAR returns a list containing the observed data (x), the true cluster membership (z), the complete data (xfull), the cluster membership given by the Baye's rule (zhat), the empirical rates of misclassification (meanerrorclass) and missngness (meanmiss).
##' @references Clustering Data with Non-Ignorable Missingness using Semi-Parametric Mixture Models, Marie Du Roy de Chaumaray and  Matthieu Marbac <arXiv:2009.07662>.
##'
##' @examples
##' set.seed(123)
##' # Data generation
##' ech <- rMNAR(n=100, K=3, d=3, delta=2, gamma=1)
##' # Head of the observed data
##' head(ech$x)
##' # Table of the cluster memberships
##' table(ech$z)
##' # Empirical rate of misclassification
##' ech$meanerrorclass
##' # Empirical rate of missingness
##' ech$meanmiss
##'
##' @export
rMNAR <- function(n, K, d=3, delta=3, gamma=1, law="gauss", linkmissing="logit-X"){
  pobs <- function(u, param) dnorm(u-param)
  xfull <- x <- matrix(NA, n, d)
  if (law == "student"){
    xfull <- x <- matrix(rt(n * d, 3), n, d)
    pobs <- function(u, param) dt(u-param, 3)
  }else if (law == "laplace"){
    xfull <- x <- matrix(rlaplace(n * d, m=0, s=1), n, d)
    pobs <- function(u, param) dlaplace(u - param,  m=0, s=1)
  }else if (law == "skewgauss"){
    xfull <- x <- matrix(rsn(n * d, xi = 0, omega = 1, alpha = 2), n, d)
    pobs <- function(u, param) dsn(u - param, xi = 0, omega = 1, alpha = 2)
  }else if (law == "gauss"){
    xfull <- x <- matrix(rnorm(n * d), n, d)
    pobs <- function(u, param) dnorm(u-param)
  }else{
    stop("law must be gauss, student, laplace or skewgauss")
  }
  missingness <- function(u, k, gamma) 1/(1 + exp(gamma))
  if (linkmissing == "MCAR"){
    missingness <- function(u, k, gamma) 1/(1 + exp(gamma))
  }else if (linkmissing == "logit-Z"){
    missingness <- function(u, k, gamma) 1/(1 + exp(gamma + 2 * k))
  }else if (linkmissing == "logit-X"){
    missingness <- function(u, k, gamma) 1/(1 + exp(gamma + u))
  }else if (linkmissing == "censoring"){
    missingness <- function(u, k, gamma) (gamma)>u
  }else{
    stop("linkmissing must be MCAR, logit-Z, logit-X or censoring")
  }
  z <- sample(1:K, n, replace=TRUE, prob =c(1/2, rep(1/(2*(K-1)), K-1)))
  for (k in unique(z)){
    who <- which(z==k)
    for (j in 1:d){
      xfull[who,  j] <- x[who, j] <- x[who, j] +  ((j-k)%%K == 0) * delta
    }
    for (j in 1:d){
      p.missing <- missingness(xfull[who,j], k, gamma)
      add.missing <- (runif(length(who)) < p.missing)
      if (any(add.missing)){
        x[who[which(add.missing==TRUE)],j] <- NA
      }
    }
  }
  logproba.class <- matrix(log(pi), n, K, byrow = TRUE)
  for (k in 1:K){
    for (j in 1:d){
      param <- 0 + ((j-k)%%K == 0) * delta
      isobs <- !is.na(x[,j])
      logproba.class[which(isobs), k] <- log(pobs(x[which(isobs), j], param)) +  log(1-missingness(x[which(isobs), j], k, gamma)) + logproba.class[which(isobs), k]
      integrand <- function(u) pobs(u, param) * missingness(u, k, gamma)
      res <- integrate(integrand, -Inf, Inf)
      logproba.class[which(isobs==FALSE), k] <-  logproba.class[which(isobs==FALSE), k] + log(res$value)
    }
  }
  zhat=apply(logproba.class, 1, which.max)
  list(x=x,
       z=z,
       xfull=xfull,
       zhat=zhat,
       meanerrorclass=1-sum(diag(table(z,zhat)))/n,
       meanmiss=mean(is.na(x)))
}
