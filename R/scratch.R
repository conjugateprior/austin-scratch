LL.assoc <- function(p, dd){
  if (is.null(p$grand))
    g <- 0 ## grand mean is optional
  else
    g <- p$grand
  
  if (is.null(p$rho))
    rho <- 1 ## rho is optional
  else
    rho <- p$rho
  
  logexpected <- g + rho * outer(p$theta, p$beta) + outer(p$alpha, p$psi, "+")
  sum(sum(dd * logexpected - exp(logexpected)))
}

stdev <- function(x){
  sqrt(var(x) * (length(x)-1) / length(x))
}

dds <- readRDS('dds.rds')
ca0 <- ca(dds)

wf0 <- austin:::wordfish(austin:::wfm(dds, word.margin=2), dir=c(1,10), control=list(tol=0.00000000001))

theta.ca0 <- ca0$rowcoord[,1]
theta.wf0 <- wf0$theta
cor(theta.ca0, theta.wf0) ## 0.999

rm <- apply(dds, 1, sum)
cm <- apply(dds, 2, sum)
eP <- rm %*% t(cm)
S <- (dds - eP)/sqrt(eP)
uvd <- svd(resid.ll0, nu=1, nv=1)

## fitted values from independence model
fitted.ll0 <- loglin(margin=c(1,2), as.matrix(dds), fit=TRUE)$fit
resid.ll0 <- (dds - fitted.ll0)/sqrt(fitted.ll0)
uvd <- svd(resid.ll0, nu=1, nv=1)
un.theta <- uvd$u
un.beta <- uvd$v * uvd$d[1]

mu.un.theta <- mean(un.theta)
rho.un.theta <- stdev(un.theta)

## unregularised wordfish model for testing purposes
wf.model <- function(dd, iterations=10){
  
  LL <- function(params){
    logexpected <- outer(params$theta, params$beta) + 
      outer(params$alpha, params$psi, "+")
    sum(sum(dd * logexpected - exp(logexpected)))
  }
  
  LL.psi.beta <- function(p, y, v.theta, v.alpha) {
    eta  <- p[1] + p[2] * v.theta + v.alpha
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  LL.alpha.theta <- function(p, y, v.psi, v.beta) {
    eta <- v.psi + v.beta * p[2] + p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
   
  LL.theta <- function(p, y, v.psi, v.beta) {
    eta <- v.psi + v.beta * p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  ## crazy start params
  N <- nrow(dd)
  V <- ncol(dd)
  params <- list(alpha=rep(0, N), 
                 psi=rep(0, V), 
                 theta=rnorm(N),
                 beta=rnorm(V))
  
  for (iter in 1:iterations){
    
    for (ii in 1:V){
      res.psi <- optim(p=c(params$psi[ii], params$beta[ii]), fn=LL.psi.beta, gr=NULL,
                       y=dd[,ii], v.theta=params$theta, v.alpha=params$alpha, 
                       method=c("BFGS"), 
                       control=list(fnscale=-1))
      if (res.psi$convergence != 0) stop(res.psi$convergence)
      params$psi[ii] <- res.psi$par[1]
      params$beta[ii] <- res.psi$par[2]
    }
    params$beta <- params$beta - mean(params$beta)
  
    res.theta <- optim(p=c(params$theta[1]), 
                       fn=LL.theta, gr=NULL,
                       y=dd[1,], v.psi=params$psi, v.beta=params$beta, 
                       method=c("BFGS"), 
                       control=list(fnscale=-1))
    if (res.theta$convergence != 0) stop(res.theta$convergence)
    params$theta[1] <- res.theta$par[1]
    
    for (ii in 2:N){
      res.alpha <- optim(p=c(params$alpha[ii], params$theta[ii]), 
                         fn=LL.alpha.theta, gr=NULL,
                         y=dd[ii,], v.psi=params$psi, v.beta=params$beta, 
                         method=c("BFGS"), 
                         control=list(fnscale=-1))
      if (res.alpha$convergence != 0) stop(res.alpha$convergence)
      params$alpha[ii] <- res.alpha$par[1]
      params$theta[ii] <- res.alpha$par[2]
    }
    params$theta <- (params$theta - mean(params$theta)) / stdev(params$theta)
    
    print(LL(params))
  }
  
  params
}

## clasest thing to ML wordscores model
ws.model <- function(dd, theta=list(), iterations=10){
  
  vals <- as.vector(unlist(theta))
  nms <- names(theta)
  fixed <- which(rownames(dd) %in% nms)
  unfixed <- which(!(rownames(dd) %in% nms))
  
  LL <- function(params){
    logexpected <- outer(params$theta, params$beta) + 
      outer(params$alpha, params$psi, "+")
    sum(sum(dd * logexpected - exp(logexpected)))
  }
  
  LL.psi.beta <- function(p, y, v.theta, v.alpha) {
    eta  <- p[1] + p[2] * v.theta + v.alpha
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  LL.alpha.theta <- function(p, y, v.psi, v.beta) {
    eta <- v.psi + v.beta * p[2] + p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  LL.theta <- function(p, y, v.psi, v.beta) {
    eta <- v.psi + v.beta * p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  LL.alpha <- function(p, y, v.psi, v.beta, fixed.theta) {
    eta <- v.psi + v.beta * fixed.theta + p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  N <- nrow(dd)
  V <- ncol(dd)

  ## hack together some starting values for theta
  theta <- rnorm(N)
  names(theta) <- rownames(dd) 
  theta[nms] <- vals
  theta <- (theta - mean(theta)) / stdev(theta)
  
  params <- list(alpha=rep(0, N), 
                 psi=rep(0, V), 
                 theta=theta,
                 beta=rnorm(V))
  
  for (iter in 1:iterations){
    
    for (ii in 1:V){
      res.psi <- optim(p=c(params$psi[ii], params$beta[ii]), fn=LL.psi.beta, gr=NULL,
                       y=dd[,ii], v.theta=params$theta, v.alpha=params$alpha, 
                       method=c("BFGS"), 
                       control=list(fnscale=-1))
      if (res.psi$convergence != 0) stop(res.psi$convergence)
      params$psi[ii] <- res.psi$par[1]
      params$beta[ii] <- res.psi$par[2]
    }
    params$beta <- params$beta - mean(params$beta)
    
    ## optimise theta and alpha (or just theta if ii==1)
    for (ii in unfixed){
      if (ii == 1){
        res.theta <- optim(p=c(params$theta[1]), 
                         fn=LL.theta, gr=NULL,
                         y=dd[1,], v.psi=params$psi, v.beta=params$beta, 
                         method=c("BFGS"), 
                         control=list(fnscale=-1))
        if (res.theta$convergence != 0) stop(res.theta$convergence)
        params$theta[1] <- res.theta$par[1]
      } else {
        res.alpha <- optim(p=c(params$alpha[ii], params$theta[ii]), 
                         fn=LL.alpha.theta, gr=NULL,
                         y=dd[ii,], v.psi=params$psi, v.beta=params$beta, 
                         method=c("BFGS"), 
                         control=list(fnscale=-1))
        if (res.alpha$convergence != 0) stop(res.alpha$convergence)
        params$alpha[ii] <- res.alpha$par[1]
        params$theta[ii] <- res.alpha$par[2]
      }
    }
    ## optimise just alphas on fixed
    for (ii in fixed){
      if (ii != 1){
        res.fixed.theta <- optim(p=c(params$alpha[ii]), 
                           fn=LL.alpha, gr=NULL,
                           y=dd[ii,], fixed.theta=params$theta[ii], 
                           v.psi=params$psi, v.beta=params$beta, 
                           method=c("BFGS"), 
                           control=list(fnscale=-1))
        if (res.fixed.theta$convergence != 0) stop(res.fixed.theta$convergence)
        params$alpha[ii] <- res.fixed.theta$par[1]
      }
    }
    ## normalise
    params$theta <- (params$theta - mean(params$theta)) / stdev(params$theta)
    
    print(LL(params))
  }
  
  ## back transform the thetas?
  
  lin <- lm(vals ~ params$theta[fixed])
  print(summary(lin))
  params$rescaled.theta <- params$theta * coef(lin)[2] + coef(lin)[1]
  
  params
}



rho.model <- function(dd, theta, beta, iterations=10){
  
  ## fixed in this model
  outer.theta.beta <- outer(theta, beta)
  
  LL <- function(params){
    logexpected <- params$rho * outer.theta.beta + 
      outer(params$alpha, params$psi, "+")
    sum(sum(dd * logexpected - exp(logexpected)))
  }
  
  LL.psi <- function(p, y, v.theta, v.alpha, s.beta, s.rho) {
    eta  <- p[1] + s.rho * s.beta * v.theta + v.alpha
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
    
  LL.alpha <- function(p, y, v.psi, v.beta, s.theta, s.rho) {
    eta <- v.psi + s.rho * v.beta * s.theta + p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  LL.rho <- function(p, Y, v.alpha, v.psi) {
    eta  <- outer(v.alpha, v.psi, "+") + p[1] * outer.theta.beta
    ll <- sum(sum(eta * Y - exp(eta)))
    return(ll)
  }
  
  ## crazy start params
  params <- list(alpha=rep(0, length(theta)), 
                 psi=rep(0, length(beta)), 
                 rho=1)
  
  ll <- LL(params)
  ll.seq <- c(ll)
  
  old.ll <- ll-1 ## guarantee one iteration
  iter <- 0
  while ((ll-old.ll) > tol){
  
  for (ii in 1:length(beta)){
    res.psi <- optim(p=c(params$psi[ii]), fn=LL.psi, gr=NULL,
      y=dd[,ii], v.theta=theta, v.alpha=params$alpha, 
      s.beta=beta[ii], s.rho=params$rho,
      method=c("BFGS"), 
      control=list(fnscale=-1))
    if (res.psi$convergence != 0) stop(res.psi$convergence)
    params$psi[ii] <- res.psi$par[1]
  }
  
  for (ii in 2:length(theta)){
    res.alpha <- optim(p=c(params$alpha[ii]), fn=LL.alpha, gr=NULL,
      y=dd[ii,], v.psi=params$psi, v.beta=beta, 
      s.theta=theta[ii], s.rho=params$rho, 
      method=c("BFGS"), 
      control=list(fnscale=-1))
    if (res.alpha$convergence != 0) stop(res.alpha$convergence)
    params$alpha[ii] <- res.alpha$par[1]
  }
  
  res.rho <- optim(p=c(params$rho), fn=LL.rho, gr=NULL,
    Y=dd, v.alpha=params$alpha, v.psi=params$psi,
    method=c("BFGS"), 
    control=list(fnscale=-1))
  if (res.rho$convergence != 0) stop(res.rho$convergence)
  params$rho <- res.rho$par[1]

  old.ll <- ll
  ll <- LL(params)
  iter <- iter + 1
  cat(paste0(iter, "\tLL=", ll, "\n"))
  ll.seq <- c(ll.seq, ll)
  }
    
  list(theta=theta, beta=beta, rho=params$rho,
       psi=params$psi, alpha=params$alpha, ll.seq=ll.seq)
}

## now a model that fixes only row and one that fixes only columns

wf2assoc <- function(mod){
  stdev <- function(x){
    sqrt(var(x) * (length(x)-1) / length(x))
  }
  
  ## interaction terms
  m.beta <- mean(mod$beta)	
  s.beta <- stdev(mod$beta)
  alpha <- mod$alpha + mod$theta * m.beta
  rho <- s.beta
  beta <- (mod$beta - m.beta)/s.beta
  ## main effect terms
  m.alpha <- mean(alpha)
  m.psi <- mean(mod$psi)
  alpha <- alpha - m.alpha
  psi <- mod$psi - m.psi
  ## constant term
  grand <- m.alpha + m.psi
  
  list(grand=grand, alpha=alpha, psi=psi, theta=mod$theta, beta=beta, rho=rho)
}

assoc2wf <- function(mod){
  alpha <- mod$alpha + mod$grand
  a1 <- alpha[1]
  alpha <- alpha - a1
  psi <- mod$psi + a1
  beta <- mod$beta * mod$rho
  
  # since WF is not quite identified we may as well leave beta centred
  list(alpha=alpha, psi=psi, theta=mod$theta, beta=beta)	
}

## unregularised wordfish model for testing purposes

alphafish.model <- function(dd, tol=0.00001){
  
  LL <- function(params){
    logexpected <- outer(params$theta, params$beta) + 
      outer(params$alpha, params$psi, "+")
    sum(sum(dd * logexpected - exp(logexpected)))
  }
  
  LL.psi.beta <- function(p, y, v.theta, v.alpha) {
    eta  <- p[1] + p[2] * v.theta + v.alpha
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
    
  LL.alpha.theta <- function(p, y, v.psi, v.beta) {
    eta <- v.psi + v.beta * p[2] + p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  LL.alpha <- function(p, y, v.psi, v.beta, fixed.theta) {
    eta <- v.psi + v.beta * p[1]
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  ## crazy start params
  N <- nrow(dd)
  V <- ncol(dd)
  params <- list(alpha=rep(0, N), 
                 psi=rep(0, V), 
                 theta=rnorm(N),
                 beta=rnorm(V))
  
  alp <- log(rowSums(dd)/sum(dd))
  psi <- log(colSums(dd))
  m.alp <- mean(alp)
  alp <- alp - m.alp
  psi <- psi + m.alp
  params$alpha <- alp
  params$psi <- psi
  
  camod <- ca(dd) ## fix me to sign flip as necessary
  params$beta <- camod$colcoord[,1] * camod$sv[1]
  params$theta <- camod$rowcoord[,1]
  
  ll <- LL(params)
  ll.seq <- c(ll)

  old.ll <- ll-1 ## guarantee one iteration
  iter <- 0
  while ((ll-old.ll) > tol){
    
    for (ii in 1:V){
      res.psi <- optim(p=c(params$psi[ii], params$beta[ii]), fn=LL.psi.beta, gr=NULL,
                       y=dd[,ii], v.theta=params$theta, v.alpha=params$alpha, 
                       method=c("BFGS"), 
                       control=list(fnscale=-1))
      if (res.psi$convergence != 0) stop(res.psi$convergence)
      params$psi[ii] <- res.psi$par[1]
      params$beta[ii] <- res.psi$par[2]
    }
    params$beta <- params$beta - mean(params$beta)
    
    for (ii in 1:N){
      res.alpha <- optim(p=c(params$alpha[ii], params$theta[ii]), 
                         fn=LL.alpha.theta, gr=NULL,
                         y=dd[ii,], v.psi=params$psi, v.beta=params$beta, 
                         method=c("BFGS"), 
                         control=list(fnscale=-1))
      if (res.alpha$convergence != 0) stop(res.alpha$convergence)
      params$alpha[ii] <- res.alpha$par[1]
      params$theta[ii] <- res.alpha$par[2]
    }
    params$alpha <- params$alpha - mean(params$alpha)
    params$theta <- (params$theta - mean(params$theta)) / stdev(params$theta)
    
    old.ll <- ll
    ll <- LL(params)
    iter <- iter + 1
    cat(paste0(iter, "\tLL=", ll, "\n"))
    ll.seq <- c(ll.seq, ll)
  }
  params$ll.seq <- ll.seq
  
  params
}
