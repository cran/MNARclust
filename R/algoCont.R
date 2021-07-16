algo_NP_smoothPDFCont <- function(ech, param){
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
    for (j in ech$contvbles){
        who <- which(ech$r[,j])
        se <- seq(min(na.omit(ech$xobs[,j])) - 5 * param$band[j], max(na.omit(ech$xobs[,j]))+ 5 * param$band[j], length.out = 1000)
        size <- (max(se) - min(se)) / length(se)
        out[who,k] <- out[who,k] + as.numeric(obj3Cpp(se, xj=ech$xobs[who,j], weightsR=param$weights[who,k],band=param$band[j])) * size
    }
  }
  out
}

algo_NP_MstepCont <- function(ech, param){
  n <- nrow(param$weights)
  param$pi <- colSums(param$weights) / ech$n
  for (k in 1:ncol(param$weights)){
      param$rho[k,] <-colSums(ech$r * param$weights[,k]) / sum(param$weights[,k])
  }
  param
}

algo_NP_onealgoCont <- function(ech, param, tol){
  logSmoothPDFwithZ <- algo_NP_smoothPDFCont(ech, param)
  logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
  logSmoothlike <- sum(logSmoothPDF)
  repeat{
    param$weights <- exp(sweep(logSmoothPDFwithZ, 1, logSmoothPDF, "-"))
    param <- algo_NP_MstepCont(ech, param)
    logSmoothPDFwithZ <- algo_NP_smoothPDFCont(ech, param)
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
