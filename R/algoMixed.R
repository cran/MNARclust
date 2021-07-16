algo_NP_smoothPDFMixed <- function(ech, param){
  K <- length(param$pi)
  out <- matrix(log(param$pi), nrow(ech$r), K, byrow = TRUE)
  for (k in 1:K){
    for (j in ech$contvbles){
      who <- which(ech$r[,j])
      out[who,k] <- out[who,k] + log(param$rho[k,j])
      if (any(ech$r[,j]==0)){
        out[which(ech$r[,j]==0),k] <-out[which(ech$r[,j]==0),k] + log(1-param$rho[k,j])
      }
      se <- seq(min(na.omit(ech$xobs[,j])) - 5 * param$band[j], max(na.omit(ech$xobs[,j]))+ 5 * param$band[j], length.out = 1000)
      size <- (max(se) - min(se)) / length(se)
      out[who,k] <- out[who,k] + as.numeric(obj3Cpp(se, xj=ech$xobs[who,j], weightsR=param$weights[who,k],band=param$band[j])) * size
    }
  }
  for (k in 1:K){
    for (j in ech$catvbles){
      out[ech$r[,j],k] <- out[ech$r[,j],k] + log(param$rho[k,j])
      if (any(ech$r[,j]==0)){
        out[which(ech$r[,j]==0),k] <-out[which(ech$r[,j]==0),k] + log(1-param$rho[k,j])
      }
      for (h in 1:nlevels(ech$xobs[,j])){
        who <- which(ech$xobs[,j]==levels(ech$xobs[,j])[h])
        out[who, k] <- out[who, k]  + log(param$alpha[[j]][k,h])
      }
    }
  }
  out
}

algo_NP_MstepMixed <- function(ech, param){
  n <- nrow(param$weights)
  param$pi <- colSums(param$weights) / ech$n
  for (k in which(param$pi!=0)){
    param$rho[k,] <- colSums(ech$r * param$weights[,k]) / sum(param$weights[,k])
  }
  for (j in ech$catvbles){
    for (h in 1:nlevels(ech$xobs[,j])){
      who <- which(ech$xobs[,j]==levels(ech$xobs[,j])[h])
      param$alpha[[j]][,h] <-  colSums(param$weights[who,,drop=FALSE])
    }
    for (k in 1:nrow(param$alpha[[j]])){
      if (sum(param$alpha[[j]][k,])>0){
        param$alpha[[j]][k,] <- param$alpha[[j]][k,] / sum(param$alpha[[j]][k,])
      }else{
        param$alpha[[j]][k,] <- 1/ncol(param$alpha[[j]])
      }
    }
  }
  param
}

algo_NP_onealgoMixed <- function(ech, param, tol){
  logSmoothPDFwithZ <- algo_NP_smoothPDFMixed(ech, param)
  logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
  logSmoothlike <- sum(logSmoothPDF)
  repeat{
    param$weights <- exp(sweep(logSmoothPDFwithZ, 1, logSmoothPDF, "-"))
    param <- algo_NP_MstepMixed(ech, param)
    logSmoothPDFwithZ <- algo_NP_smoothPDFMixed(ech, param)
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
