algo_NP_smoothPDFCat <- function(ech, param){
  K <- length(param$pi)
  out <- matrix(log(param$pi), nrow(ech$r), K, byrow = TRUE)
  for (k in 1:K){
    valmissing <- rep(0, ncol(param$rho))
    if (any(param$rho[k,]<1)){
      where <- which(param$rho[k,]<1)
      valmissing[where] <- log(1 - param$rho[k,where])
    }
    out[,k] <- out[,k] + rowSums(sweep(x = ech$r, MARGIN = 2, STATS = log(param$rho[k,]), FUN =  "*") +
                                   sweep(x = 1 - ech$r, MARGIN = 2, STATS = valmissing, FUN = "*"))
    for (j in ech$catvbles){
      for (h in 1:nlevels(ech$xobs[,j])){
        who <- which(ech$xobs[,j]==levels(ech$xobs[,j])[h])
        out[who, k] <- out[who, k]  + log(param$alpha[[j]][k,h])
      }

    }
  }
  out
}

algo_NP_MstepCat <- function(ech, param){
  n <- nrow(param$weights)
  param$pi <- colSums(param$weights) / ech$n
  for (k in 1:ncol(param$weights)){
    param$rho[k,] <-colSums(ech$r * param$weights[,k]) / sum(param$weights[,k])
  }
  for (j in ech$catvbles){
    for (h in 1:nlevels(ech$xobs[,j])){
      who <- which(ech$xobs[,j]==levels(ech$xobs[,j])[h])
      param$alpha[[j]][,h] <-  colSums(param$weights[who,,drop=FALSE])
    }
    for (k in 1:nrow(param$alpha[[j]])){
      param$alpha[[j]][k,] <- param$alpha[[j]][k,] / sum(param$alpha[[j]][k,])
    }
  }
  param
}

algo_NP_onealgoCat <- function(ech, param, tol){
  logSmoothPDFwithZ <- algo_NP_smoothPDFCat(ech, param)
  logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
  logSmoothlike <- sum(logSmoothPDF)
  repeat{
    param$weights <- exp(sweep(logSmoothPDFwithZ, 1, logSmoothPDF, "-"))
    param <- algo_NP_MstepCat(ech, param)
    logSmoothPDFwithZ <- algo_NP_smoothPDFCat(ech, param)
    logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
    prec <- logSmoothlike
    logSmoothlike <- sum(logSmoothPDF)
    if (is.nan(logSmoothlike)){
      logSmoothlike <- prec <- -Inf
      break
    }
    if ( (logSmoothlike - prec) < tol){
      break
    }
  }
  list(param = param, logSmoothlike = logSmoothlike, completedSmoothloglike = sum(apply(logSmoothPDFwithZ, 1, max)), logSmoothPDFwithZ = logSmoothPDFwithZ, prec = prec, zhat = apply(param$weights, 1, which.max))
}


