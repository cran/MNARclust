builddata <- function(x){
  contvbles <- catvbles <- NULL
  if (is.data.frame(x)){
    contvbles <- which(sapply(x, is.numeric))
    catvbles <- which(!sapply(x, is.numeric))
  }else{
    if (is.numeric(x[,1]))
      contvbles <- 1:ncol(x)
    else
      catvbles <- 1:ncol(x)
  }
  list(xobs = x, r = !is.na(x), contvbles=contvbles, catvbles=catvbles, d=ncol(x), n=nrow(x))
}


logRowSums <- function(u){
  tmp <- apply(u, 1, max)
  tmp + log(rowSums(exp(sweep(u, 1, tmp, "-"))))
}

start.paramCont <- function(ech, K, band){
  centers <- sample(1:ech$n, K)
  pi <- rep(1/K, K)
  rho <- matrix(runif(K * ech$d)/3 + .33, K, ech$d)
  log.weights <- matrix(log(pi), ech$n, K, byrow = TRUE)
  for (k in 1:K){
    log.weights <- log.weights + rowSums(sweep(x = ech$r, MARGIN = 2, STATS = log(rho[k,]), FUN =  "*") + sweep(x = 1 - ech$r, MARGIN = 2, STATS = log(1 - rho[k,]), FUN = "*"))
    for (j in 1:ech$d){
      mu <- mean(ech$xobs[,j], na.rm = TRUE)
      sd <- sd(ech$xobs[,j], na.rm = TRUE)
      if (ech$r[centers[k], j])  mu <- ech$xobs[centers[k],j]
      if (sum(ech$r[,j])>0){
        who.tmp <- which(ech$r[,j]==1)
        log.weights[who.tmp,k] <- log.weights[who.tmp,k] + dnorm(x = ech$x[who.tmp,j], mean = mu, sd = sd, log = TRUE)
      }
    }
  }
  weights <- exp(sweep(log.weights, 1, logRowSums(log.weights), "-"))
  list(pi = pi, rho = rho, weights = weights, band = band)
}

start.paramMixed <- function(ech, K, band){
  centers <- sample(1:ech$n, K)
  pi <- rep(1/K, K)
  rho <- matrix(runif(K * ech$d)/3 + .33, K, ech$d)
  log.weights <- matrix(log(pi), ech$n, K, byrow = TRUE)
  alpha <- list()
  length(alpha) <- ech$d
  for (j in ech$catvbles){
    alpha[[j]] <- matrix(NA, K, nlevels(ech$xobs[,j]))
  }
  for (k in 1:K){
    log.weights <- log.weights + rowSums(sweep(x = ech$r, MARGIN = 2, STATS = log(rho[k,]), FUN =  "*") + sweep(x = 1 - ech$r, MARGIN = 2, STATS = log(1 - rho[k,]), FUN = "*"))
    for (j in ech$contvbles){
      mu <- mean(ech$xobs[,j], na.rm = TRUE)
      sd <- sd(ech$xobs[,j], na.rm = TRUE)
      if (ech$r[centers[k], j])  mu <- ech$xobs[centers[k],j]
      if (sum(ech$r[,j])>0){
        who.tmp <- which(ech$r[,j]==1)
        log.weights[who.tmp,k] <- log.weights[who.tmp,k] + dnorm(x = ech$x[who.tmp,j], mean = mu, sd = sd, log = TRUE)
      }
    }
    for (j in ech$catvbles){
      tmpalpha <- runif(nlevels(ech$xobs[,j]))
      if (ech$r[centers[k], j]) tmpalpha[ech$xobs[centers[k],j]] <-  tmpalpha[ech$xobs[centers[k],j]] + runif(1)
      tmpalpha <-  tmpalpha / sum(tmpalpha)
      alpha[[j]][k,] <- tmpalpha
      if (sum(ech$r[,j])>0){
        who.tmp <- which(ech$r[,j]==1)
        log.weights[who.tmp,k] <- log.weights[who.tmp,k] + log(tmpalpha[ech$xobs[,j]])[who.tmp]
      }
    }
  }
  weights <- exp(sweep(log.weights, 1, logRowSums(log.weights), "-"))
  list(pi = pi, rho = rho, weights = weights, band = band, alpha = alpha)
}

start.paramCat <- function(ech, K, band){
  centers <- sample(1:ech$n, K)
  pi <- rep(1/K, K)
  rho <- matrix(runif(K * ech$d)/3 + .33, K, ech$d)
  log.weights <- matrix(log(pi), ech$n, K, byrow = TRUE)
  alpha <- list()
  length(alpha) <- ech$d
  for (j in ech$catvbles){
    alpha[[j]] <- matrix(NA, K, nlevels(ech$xobs[,j]))
  }

  for (k in 1:K){
    log.weights <- log.weights + rowSums(sweep(x = ech$r, MARGIN = 2, STATS = log(rho[k,]), FUN =  "*") + sweep(x = 1 - ech$r, MARGIN = 2, STATS = log(1 - rho[k,]), FUN = "*"))
    for (j in ech$catvbles){
      tmpalpha <- runif(nlevels(ech$xobs[,j]))
      if (ech$r[centers[k], j]) tmpalpha[ech$xobs[centers[k],j]] <-  tmpalpha[ech$xobs[centers[k],j]] + runif(1)
      tmpalpha <-  tmpalpha / sum(tmpalpha)
      alpha[[j]][k,] <- tmpalpha
      if (sum(ech$r[,j])>0){
        who.tmp <- which(ech$r[,j]==1)
        log.weights[who.tmp,k] <- log.weights[who.tmp,k] + log(tmpalpha[ech$xobs[,j]])[who.tmp]
      }
    }
  }
  weights <- exp(sweep(log.weights, 1, logRowSums(log.weights), "-"))
  list(pi = pi, rho = rho, weights = weights, band = band, alpha = alpha)
}

band.default <- function(x){
  out <- rep(nrow(x)**(-1/5), ncol(x))
  for (j in 1:ncol(x)){
    if (is.numeric(x[,j])) out[j] <- sd(x[,j], na.rm=T) * out[j]
  }
  out
}


