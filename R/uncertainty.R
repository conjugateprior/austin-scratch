##' Parametric bootstrap for a wordfish model
##'
##' Runs a parametric bootstrap on a wordfish model, i.e. asuming that the model is
##' true.  This is asymptotically equivalent to the Hessian-derived standard errors
##' provided by wordfish.  The function creates L new word frequency matrices by
##' sampling from the fitted model and returns a matrix of estimated thetas (row scores)
##' To summarise these as standard errors take the standard deviation
##' or compute quantiles of the returned matrix across rows.
##'
##' The default value of L is unrealistically small, particularly if quantile
##' statistics are used to summarise the bootstrap samples.
##'
##' Note: With some probability words that occur in the original word frequency matrix
##' will not appear in any bootstrapped word frequency matrix.  These are quietly
##' discarded before fitting a model to the bootstrapped counts.
##' 
##' @title Parametric bootstrap for a wordfish model 
##' @param m wordfish model
##' @param L number of boostrap samples  
##' @param fit.fun engine used to fit the model.  Should be the same as the engine used to fit m
##' @param verbose whether to generate a running commentary
##' @return For D documents, a D by L matrix of bootstrapped theta values (row scores)  
##' @export
##' @author Will Lowe
parametric.bootstrap.se <- function(m, L=50, fit.fun='classic.wordfish', verbose=FALSE){
  D <- length(m$theta)
  V <- length(m$beta)
  mtheta <- matrix(0, nrow=D, ncol=L)
  logexpected <- outer(m$theta, m$beta) + outer(m$alpha, m$psi, "+")
  lambda <- exp(logexpected)
  for (l in 1:L){
    if (verbose) cat(paste("iter:", l, "\n"))
    mat <- matrix(rpois(rep(1, D*V), as.vector(lambda)), nrow=D)
    rownames(mat) <- m$docs
    colnames(mat) <- m$words
    newparams <- wordfish(mat, dir=m$dir, tol=m$tol, sigma=m$sigma, params=m,
                          fit.fun=fit.fun, verbose=verbose)
    mtheta[,l] <- newparams$theta
  }  
  return(mtheta)
}

