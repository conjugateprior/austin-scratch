##' Tidy up count matrix before estimation
##'
##' This function removes words that do not occur in the
##' corpus.  Mostly useful for when bootstrap functions
##' make degenerate word count matrices.
##' 
##' @title Tidy up count matrix before estimation 
##' @param wfm word count matrix: document as rows, words as columns
##' @param verbose whether to offer running commentary
##' @return a word count matrix
##' @export
##' @author Will Lowe
prepare.data <- function(wfm, verbose){
  ## strip out any words that for whatever reason do not occur in any document
  origV <- nrow(wfm)
  uninformative <- apply(wfm, 2, sum) == 0
  if (sum(uninformative)>0 && verbose)
    cat('Ejecting words that never occur: ', paste((1:origV)[uninformative], sep=','), '\n')
  good.words <- colnames(wfm)[!uninformative] ## not indices, words
  return(wfm[,good.words])
}

##' Fit a wordfish model S&P style
##'
##' This engine fits the Slapin and Proksch parameterisation of a regularized version
##' of Goodman's RC model by coordinate ascent,
##' alternating the fitting of row scores theta (with main effects alpha)
##' and column scores beta (with main effects psi).
##' This version of the model is identified by imposing 1) a zero mean
##' and unit variance constraint on the row scores (theta), 2) a fixed zero for the first
##' row effect (alpha), and 3) a directional constraint on row scores:
##' theta[dir[1]] < theta[dir[2]].
##' 
##' There are some minor differences to previous implemenentations:
##' Unlike Goodman's RC model, the beta is the MAP estimated computed with a
##' mean zero, standard deviation sigma Normally distrubuted prior.
##' Unlike Slapin and Proksch's implementation, 'standard errors' for theta are computed from the
##' Hessian, assuming other parameters fixed, rather than using a parametric bootstrap.
##' This is quicker, asymptotically equivalent, and tends to give very similar results.
##' 
##' @title Fit a wordfish model in S&P parameterisation
##' @param wfm a word count matrix: documents are rows and columns are words
##' @param dir identification constraint: theta[dir[1]] will be less than theta[dir[2]]
##' @param tol the log likelihood amount less than which parameter estimation will stop
##' @param sigma ridge regression regularization parameter
##' @param params an initialized set of parameter to start estimation from
##' @param verbose whether to give running commentary
##' @export
##' @return a fitted wordfish model
##' @author Will Lowe
classic.wordfish <- function(wfm, dir, tol, sigma, params, verbose=FALSE){

  LL.reg <- function(params, wfm){
    logexpected <- outer(params$theta, params$beta) +
      outer(params$alpha, params$psi, "+")
    sum(sum(wfm * logexpected - exp(logexpected))) - sum(0.5*(params$beta^2/sigma^2))
  }
  
  ## estimate all psi and beta
  LL.psi.beta <- function(p, y, theta, alpha, sigma) {
    beta <- p[1]
    psi  <- p[2]
    eta  <- psi + beta*theta + alpha
    ll <- sum(eta*y - exp(eta)) - 0.5*(beta^2/sigma^2)
    return(ll)
  }
  
  DLL.psi.beta <- function(p, y, theta, alpha, sigma){
    beta <- p[1]
    psi <- p[2]
    mu <- exp(psi + beta*theta + alpha)
    dll <- c(sum(theta * (y-mu)) - beta/sigma^2,
             sum(y-mu))
    return(dll)
  }
  
  ## estimate theta[1]
  LL.first.theta <- function(p, y, beta, psi) {
    theta <- p[1]
    eta <- psi + beta*theta   # alpha=0
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  LL.alpha.theta <- function(p, y, beta, psi) {
    theta <- p[1]
    alpha <- p[2]
    eta <- psi + beta*theta + alpha
    ll <- sum(eta*y - exp(eta))
    return(ll)
  }
  
  DLL.alpha.theta <- function(p, y, beta, psi){
    theta <- p[1]
    alpha <- p[2]
    mu <- exp(psi + beta*theta + alpha)
    dll <- c(sum(beta * (y-mu)),
             sum(y-mu))
    return(dll)
  }
  
  D  <- nrow(wfm)
  V  <- ncol(wfm)  
  oldll <- -Inf
  ll <- LL.reg(params, wfm)
  ll.path <- ll ## keep track
  iter <- 0
  if (verbose) cat(iter, "\tLL(reg):", ll, "\n")

  ## check that we're not going backwards but not yet done
  prog <- function(oldll, newll, tol){
    (newll-oldll) > 0 && (newll-oldll) > tol
  }

  iter <- 1
  while ( prog(oldll,ll,tol) ) {
    if (verbose) cat(iter, "\t")
    
    resa <- optim(p=c(params$theta[1]), fn=LL.first.theta, gr=NULL,
                  y=as.numeric(wfm[1,]), beta=params$beta, psi=params$psi,
                  method=c("BFGS"), control=list(fnscale=-1))
    params$theta[1] <- resa$par[1]
    params$alpha[1] <- 0  ## just in case
    if (resa$convergence != 0)
      warning("optim failed while estimating theta[ 1 ]:", colnames(wfm)[1])
    
    for (i in 2:D) {
      resa <- optim(par=c(params$theta[i], params$alpha[i]),
                    fn=LL.alpha.theta, gr=DLL.alpha.theta,
                    y=as.numeric(wfm[i,]), beta=params$beta, psi=params$psi,
                    method=c('BFGS'), control=list(fnscale=-1))
      params$theta[i] <- resa$par[1]
      params$alpha[i] <- resa$par[2]
      if (resa$convergence!=0)
        warning("optim failed while estimating theta[", i, "] and alpha[", i, "]:", colnames(wfm)[i])
    }

    m.theta <- mean(params$theta)
    s.theta <- sqrt(mean((params$theta-m.theta)**2))
    params$theta <- (params$theta - m.theta)/s.theta

    ## oldbeta <- params$beta; params$beta  <- params$beta * sd(params$theta);
    ## params$psi <- params$psi + oldbeta*thetabar    
    if (params$theta[dir[1]] > params$theta[dir[2]])
      params$theta <- params$theta*(-1)
    
    for (j in 1:V) {
      resa <- optim(par=c(params$beta[j], params$psi[j]), fn=LL.psi.beta, gr=DLL.psi.beta,
                    y=wfm[,j], theta=params$theta, alpha=params$alpha, sigma=sigma,
                    method=c('BFGS'), control=list(fnscale=-1)
                    )
      params$beta[j] <- resa$par[1]
      params$psi[j] <- resa$par[2]
      if (resa$convergence!=0)
        warning("optim failed while estimating beta[", j, "] and psi[", j, "]")
    }

    new.ll <- LL.reg(params, wfm)
    if (verbose) {
      cat("LL(reg):", ll, "\t(diff:", (new.ll-ll), ")\n")
      flush.console()
    }
    oldll <- ll
    ll <- new.ll
    ll.path <- c(ll.path, ll)
    iter <- iter + 1
  }
  if (verbose)
    cat("Log likelihood:", (LL.reg(params, wfm) + sum(0.5*(params$beta^2/sigma^2))), "\n")
  
  model <- list(dir=dir,theta=params$theta, theta=params$theta, alpha=params$alpha,
                beta=params$beta, psi=params$psi, docs=rownames(wfm),
                words=colnames(wfm), sigma=sigma,
                ll=ll, ll.path=ll.path, data=wfm)

  return(model)  ## just a list of stuff at this point
}

##' Fit a wordfish model
##'
##' This engine fits the Slapin and Proksch parameterisation of a regularized version
##' of Goodman's RC model by coordinate ascent,
##' alternating the fitting of row scores theta (with main effects alpha)
##' and column scores beta (with main effects psi).
##' This version of the model is identified by imposing 1) zero mean
##' and unit variance constrants on row scores (theta), 2) a fixed zero for the first
##' row effect (alpha), and 3) a directional constraint on row scores:
##' theta[dir[1]] < theta[dir[2]].
##' 
##' There are some minor differences to previous implemenentations:
##' Unlike Goodman's RC model, beta is ridge regularised with standard deviation sigma,
##' there is no grand mean and no separate correlation parameter.
##' Unlike Slapin and Proksch's implementation, 'standard errors' for theta are computed from the
##' Hessian, assuming other parameters fixed rather than with a parametric bootstrap.
##' 
##' @title Fit a wordfish model
##' @param wfm a word count matrix, documents are rows and columns are words
##' @param dir identification constraint: theta[dir[1]] will be less than theta[dir[2]]
##' @param tol the log likelihood amount less than which parameter estimation will stop
##' @param sigma is ridge regression regularization parameter
##' @param params is an initialized set of parameter to start estimation from or NULL, which will trigger the application of init.fun
##' @param init.fun function to generate an initial set of parameters. By default: loglinear.initialize.wordfish
##' @param fit.fun the engine to use to fit the model. By default: classic.wordfish
##' @param verbose whether to give running commentary on estimation
##' @export
##' @return a fitted wordfish model
##' @author Will Lowe
wordfish <- function(wfm, dir=c(1, 10), tol=1e-06, sigma=3, params=NULL, init.fun='loglinear.initialize.wordfish', fit.fun='classic.wordfish', verbose=FALSE){

  thecall <- match.call()
  wfm <- prepare.data(wfm, verbose)

  D  <- nrow(wfm)
  V  <- ncol(wfm)
  if (is.null(params))
    params <- do.call(init.fun, list(wfm), quote=TRUE) 
  params <- do.call(fit.fun, list(wfm, dir, tol, sigma, params, verbose), quote=TRUE)

  ## SEs: switch to multinomial parameterisation to remove alpha
  new.beta <- params$beta[2:V] - params$beta[1]
  new.psi <- params$psi[2:V] - params$psi[1]
  theta.se <- rep(NA, D)
  multinomial.ll <- function(x, y){
    dmultinom(y, prob=exp(c(0, x*new.beta + new.psi)), log=TRUE)  
  }
  for (i in 1:D){
    neghess <- -hessian(multinomial.ll, params$theta[i], y=wfm[i,])
    theta.se[i] <- sqrt(1/neghess)
  }
  params$se.theta <- theta.se
  params$call <- thecall
  
  class(params) <- c("wordfish", class(params))
  return(params)
}


##' Initialize wordfish parameters
##'
##' This is mostly the Slapin and Proksch initialization code, which is
##' in turn mostly the NOMINATE code involving an SVD of a double-centred
##' word count matrix.  Use loglinear.initialize.wordfish instead.  It'll
##' be much quicker.
##' 
##' @title Initialize wordfish parameters
##' @param wfm word count matrix with documents as rows and words as columns 
##' @export
##' @return a set of wordfish parameters: alpha, theta, psi and beta
##' @author Will Lowe (from code by J. Slapin and S.-O. Proksch)
initialize.wordfish <- function(wfm){
  D             <- nrow(wfm)
  V             <- ncol(wfm)
  numword       <- rep(1:V, each=D)
  numdoc        <- rep(1:D, V)
  
  dat           <- matrix(1, nrow=V*D, ncol=3)
  dat[,1]       <- as.vector(as.matrix(wfm))
  dat[,2]       <- as.vector(numword)
  dat[,3]       <- as.vector(numdoc)
  dat           <- data.frame(dat)
  
  colnames(dat) <- c("y", "word", "doc")
  dat$word      <- factor(dat$word)
  dat$doc       <- factor(dat$doc)
  
  b0            <- log(colMeans(wfm))
  
  rmeans        <- rowMeans(wfm)
  alpha         <- log(rmeans/rmeans[1])
  
  ystar         <- log(dat$y+0.1) - alpha[dat$doc] - b0[dat$word]
  ystarm        <- matrix(ystar, D, V, byrow=FALSE)
  res           <- svd(ystarm, nu=1)
  b             <- as.numeric(res$v[,1]*res$d[1])
  
  ## normalize
  th            <- as.vector(res$u)-res$u[1,1]
  theta         <- th/sd(th)
  b1            <- b*sd(th)
  params        <- list(alpha=as.numeric(alpha), psi=as.numeric(b0), beta=b1, theta=theta)
  return(params)
}

##' Initialise wordfish parameters from a loglinear model
##'
##' This function uses a loglinear model as the basis for parameter estimates.
##' First, we use loglin to fit a saturated model, adding a half to all counts to
##' avoid numerical difficulties.  The log linear parameters
##' consist of an intercept, the marginal effects 'docs' and 'words', and a matrix of
##' interaction terms 'docs.words'.  The interaction matrix is then decomposed via SVD.
##' Starting values for theta are taken from the first column of U and beta values from the
##' first column of U (in the orientation given by R's svd function).  Starting values
##' for alpha and psi are taken from the intercept and main effect parameters.  All parameters
##' are then transformed so that theta is zero mean unit variance, alpha[1] is zero, and
##' the intercept is absorbed into the main effects as per the usual wordfish
##' parameterisation.
##' 
##' @title Initialise wordfish parameters from a loglinear model
##' @param wfm word count matrix with documents as rows and words as columns
##' @return An initial set of wordfish model parameters
##' @export
##' @author Will Lowe
loglinear.initialize.wordfish <- function(wfm){
  ## fit a cheap saturated model, mostly so we can
  ## decompose the interaction term
  mm <- loglin(as.matrix(wfm) + 0.1, margin=list(c(1,2)), param=TRUE, print=FALSE)
  alpha.s <- mm$param[[2]]
  psi.s <- mm$param[[3]]
  grand <- mm$param[[1]]
  res <- svd(mm$param[[4]])
  theta.s <- res$u[,1]
  beta.s <- res$v[,1]

  ## rejig into the wordfish parameterisation
  alpha <- alpha.s + grand
  psi <- psi.s
  alpha1 <- alpha[1]
  alpha <- alpha - alpha1 ## zero first alpha
  psi <- psi + alpha1

  m.theta <- mean(theta.s)
  sd.theta <-sqrt(mean((theta.s - m.theta)**2))
  theta <- (theta.s - m.theta)/sd.theta
  beta <- beta.s * sd.theta + m.theta

  list(alpha=alpha, psi=psi, beta=beta, theta=theta)
}

##' Orthonormalize the columns of a matrix
##'
##' This function orthogonalises the columns of a matrix u using the
##' Gram-Schmidt method.  It is essentially the orthonormalisation
##' function from the far package by Julien Damon.
##'
##' @title Orthonormalise the columns of a matrix 
##' @param u a matrix
##' @param basis whether to create a basis too
##' @param norm whether to orthonormalise too
##' @return a matrix the same size with orthogonal columns 
##' @author Will Lowe (original code by Damon Julien Guillas Serge)
orthonormalize <- function(u, basis=FALSE, norm=FALSE){
  p <- nrow(u)  # dimension of the space
  n <- ncol(u)  # number of vectors
  if (prod(abs(La.svd(u)$d)>1e-8)==0)
    stop("colinear columns")

  # if a basis is desired then u is extended to be square
  if (basis && (p > n)) {
    base <- diag(p) # the canonical basis of size p x p
    coef.proj <- crossprod(u,base) / diag(crossprod(u))
    ## calculation of base2,
    ## the orthogonal (to the span space of u) part of the base
    base2 <- base - u %*% matrix(coef.proj,nrow=n,ncol=p)
    ## norms of the vector of base2
    norm.base2 <- diag(crossprod(base2))
    base <- as.matrix(base[,order(norm.base2) > n])
    ## extention of u with the more orthogonal vector of base
    u <- cbind(u,base)
    n <- p
  }
  
  v <- u  ## start of gram-schmidt 
  if (n > 1){
    for (i in 2:n){
      coef.proj <- c(crossprod(u[,i],v[,1:(i-1)])) / diag(crossprod(v[,1:(i-1)]))
      v[,i] <- u[,i] - matrix(v[,1:(i-1)],nrow=p) %*% matrix(coef.proj,nrow=i-1)
    }
  }
  if (norm){
    coef.proj <- 1/sqrt(diag(crossprod(v)))
    v <- t(t(v) * coef.proj)
  }
  return(v)
} 

##' Fit a wordscores model
##'
##' This function computes wordscores using the algorithm
##' described in Laver and Garry (2003) and implemented as classic.wordscores.
##' To score 'virgin', i.e out-of-sample documents, use the predict function.
##' Words that do not occur in the document set are quietly removed before
##' wordscores are computed.
##' 
##' @title Fit a wordscores model
##' @param wfm a word count matrix, with documents as rows and columns as words
##' @param scores reference scores for all documents
##' @param fit.fun engine to use to fit the model
##' @param verbose whether to give a running commentary
##' @return a wordscores model
##' @export
##' @author Will Lowe
wordscores <- function(wfm, scores, fit.fun='classic.wordscores', verbose=FALSE){
  thecall <- match.call()
  wfm <- prepare.data(wfm, verbose)
  params <- do.call(fit.fun, list(wfm, scores), quote=TRUE)
  params$call <- thecall
  class(params) <- c('wordscores',  class(params)) 
  return(params)
}
##' Fit a wordscores model
##'
##' Fit a Wordscores model using the algorithm in Laver and Garry.
##' This means tossing out the words that never occur and
##' constructing wordscores.
##' Consider using correspondence analysis instead!
##' 
##' @title Fit a wordscores model 
##' @param wfm a word count matrix with documents as rows and words as columns
##' @param scores reference document scores
##' @export
##' @return a fitted wordscores model 
##' @author Will Lowe
classic.wordscores <- function(wfm, scores){
  C <- wfm[,colSums(wfm)>0]  ## just the words that occur
  F <- C / outer(rep(1, nrow(C)), colSums(C))
  pi <- as.numeric((scores %*% F) / colSums(F))
  names(scores) <- rownames(wfm) ## just in case
  names(pi) <- colnames(C)
  list(pi=pi, theta=scores, data=wfm)
}


rcm.wordfish <- function(wfm, dir, tol, sigma, params, verbose=FALSE){

  LL.reg <- function(params, wfm){
    logexpected <- params$phi * outer(params$theta, params$beta) +
      outer(params$alpha, params$psi, "+") + params$grand
    sum(sum(wfm * logexpected - exp(logexpected)))
  }
  
  LL.psi.beta <- function(params, grand, phi, theta, alpha){    
    logexpected <- phi * outer(theta, params[1:V]) +
      outer(alpha, params[(V+1):(2*V)], "+") + grand
    sum(sum(wfm * logexpected - exp(logexpected)))
  }

  LL.alpha.theta <- function(params, grand, phi, beta, psi){
    logexpected <- phi * outer(params[1:D], beta) +
      outer(params[(D+1):(2*D)], psi, "+") + grand
    sum(sum(wfm * logexpected - exp(logexpected)))
  }

  LL.grand <- function(params, phi, beta, psi, theta, alpha){
    logexpected <- phi * outer(theta, beta) +
      outer(alpha, psi, "+") + params[1]
    sum(sum(wfm * logexpected - exp(logexpected)))
  }

  LL.phi <- function(params, grand, beta, psi, theta, alpha){
    logexpected <- params[1] * outer(theta, beta) +
      outer(alpha, psi, "+") + grand
    sum(sum(wfm * logexpected - exp(logexpected)))
  }
  
  D  <- nrow(wfm)
  V  <- ncol(wfm)  
  oldll <- -Inf
  ll <- LL.reg(params, wfm)
  ll.path <- ll ## keep track
  iter <- 0
  if (verbose) cat(iter, "\tLL(reg):", ll, "\n")

  ## check that we're not going backwards but not yet done
  prog <- function(oldll, newll, tol){
    (newll-oldll) > 0 && (newll-oldll) > tol
  }

  iter <- 1
  while ( prog(oldll,ll,tol) ) {
    if (verbose) cat(iter, "\t")
        
    resa <- optim(par=c(params$theta, params$alpha),
                  fn=LL.alpha.theta, gr=NULL,
                  beta=params$beta, psi=params$psi, grand=params$grand, phi=params$phi,
                  method=c('BFGS'), control=list(fnscale=-1))
    if (resa$convergence!=0)
      warning("optim failed in alpha-theta")
    params$theta <- resa$par[1:D]
    params$alpha <- resa$par[(D+1):(2*D)]

    m.theta <- mean(params$theta)
    s.theta <- sqrt(mean((params$theta-m.theta)**2))
    params$theta <- (params$theta - m.theta)/s.theta

	m.alpha <- mean(params$alpha)
	params$alpha <- params$alpha - m.alpha

    if (params$theta[dir[1]] > params$theta[dir[2]])
      params$theta <- params$theta*(-1)
    
    resa <- optim(par=c(params$beta, params$psi), fn=LL.psi.beta, gr=NULL,
                  grand=params$grand, phi=params$phi, theta=params$theta, alpha=params$alpha,
                  method=c('BFGS'), control=list(fnscale=-1)
                  )
    if (resa$convergence!=0)
      warning("optim failed in psi-beta")
    params$beta <- resa$par[1:V]
    params$psi <- resa$par[(V+1):(2*V)]

    m.beta <- mean(params$beta)
    s.beta <- sqrt(mean((params$beta-m.beta)**2))
    params$beta <- (params$beta - m.beta)/s.beta

	m.psi <- mean(params$psi)
	params$psi <- params$psi - m.psi

    resa <- optim(par=c(params$grand), fn=LL.grand, gr=NULL,
                  beta=params$beta, psi=params$psi, phi=params$phi, theta=params$theta, alpha=params$alpha,
                  method=c('BFGS'), control=list(fnscale=-1)
                  )
    if (resa$convergence!=0)
      warning("optim failed in grand")
    params$grand <- resa$par[1]

    resa <- optim(par=c(params$phi), fn=LL.phi, gr=NULL,
                  beta=params$beta, psi=params$psi, grand=params$grand, theta=params$theta, alpha=params$alpha,
                  method=c('BFGS'), control=list(fnscale=-1)
                  )
    if (resa$convergence!=0)
      warning("optim failed in phi")
    params$phi <- resa$par[1]

    new.ll <- LL.reg(params, wfm)
    if (verbose) {
      cat("LL(reg):", ll, "\t(diff:", (new.ll-ll), ")\n")
      flush.console()
    }
    oldll <- ll
    ll <- new.ll
    ll.path <- c(ll.path, ll)
    iter <- iter + 1
  }
  if (verbose)
    cat("Log likelihood:", LL.reg(params, wfm), "\n")
  
  model <- list(dir=dir,theta=params$theta, alpha=params$alpha,
                beta=params$beta, psi=params$psi, phi=params$phi, grand=params$grand, 
                docs=rownames(wfm),
                words=colnames(wfm),
                ll=ll, ll.path=ll.path, data=wfm)

  return(model)  ## just a list of stuff at this point
}

rcm.init <- function(wfm){
	mm <- loglin(dd$Y+0.1, margin=list(c(1,2)), param=TRUE)	
	params <- list(grand=mm$param$`(Intercept)`,
		psi=mm$param$words, alpha=mm$param$docs)
	inter <- svd(mm$param$docs.words)
	params$theta <- inter$u[,1]
	params$beta <- inter$v[,1]
	params$phi <- 0 # inter$d[1]
	
	m.beta <- mean(params$beta)
    s.beta <- sqrt(mean((params$beta-m.beta)**2))
    params$beta <- (params$beta - m.beta)/s.beta
    
    m.theta <- mean(params$theta)
    s.theta <- sqrt(mean((params$theta-m.theta)**2))
    params$theta <- (params$theta - m.theta)/s.theta

	return(params)
}
	

