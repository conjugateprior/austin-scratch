##' Show a wordfish model
##'
##' Shows the fitted document positions and intervals for a fitted wordfish model.
##' @title Show a wordfish model
##' @param x wordfish model
##' @param digits how many decimal places to use
##' @param ... Not used.
##' @return Nothing. Used for printing side-effect
##' @export
##' @method print wordfish
##' @author Will Lowe
print.wordfish <- function(x, digits=max(3,getOption('digits')-3), ...){
  cat("Call:\n\t")
  print(x$call)
  cat("\nDocument Positions:\n")
  pp <- predict(x, se.fit=TRUE, interval='confidence')
  colnames(pp) <- c("Estimate", "Std. Error", "Lower", "Upper")
  print(pp, digits=digits)
}

##' Show wordfish model parameters
##'
##' Shows the parameters of a wordfish model in various equivalent parameterisations.
##' Currently these are 'poisson', the original Wordfish and RC model version, and
##' 'multinomial', the equivalent model when row totals (document lengths) are
##' assumed to be fixed.
##'
##' Note: This function does not return estimated document positions (row scores) theta. Use
##' summary for that.
##' 
##' @title Show wordfish model parameters
##' @param object fitted wordfish model
##' @param form which parameterization to return
##' @param ... Not used.
##' @return a list containing two sets of word parameters and document main effects if form='poisson'
##' @export
##' @method coef wordfish
##' @author Will Lowe
coef.wordfish <- function(object, form=c('poisson', 'multinomial'), ...){
  m <- object
  fm <- match.arg(form)
  if (fm == 'multinomial'){
    V <- length(m$beta)
    word.params <- data.frame(beta=(m$beta[2:V] - m$beta[1]),
                              psi=(m$psi[2:V] - m$psi[1]))
    rownames(word.params) <- paste(m$words[2:V], "/", m$words[1], sep='')
    d <- list(words=word.params)
  } else {
    word.params <- data.frame(beta=m$beta, psi=m$psi)
    rownames(word.params) <- m$words
    docs <- data.frame(alpha=m$alpha)
    rownames(docs) <- m$docs
    d <- list(words=word.params, docs=docs)
  }
  
  class(d) <- c('coef.wordfish', class(d))
  return(d)
}

##' Print model parameters
##'
##' Prints out the model parameters in whatever parameterisation they were asked for in.
##' @title Print model parameters
##' @param x wordfish model parameters from coef
##' @param digits how many decinmal places to show
##' @param ... Not used
##' @return Nothing. Used for printing side-effect
##' @export
##' @method print coef.wordfish
##' @author Will Lowe
print.coef.wordfish <- function(x, digits=max(3,getOption('digits')-3), ...){
  print(x$words, digits=digits)
  if (!is.null(x$docs))
    print(x$docs, digits=digits)
}

##' Plot model parameters
##'
##' Plots word parameter extracted using coef with form='poisson', in the familiar
##' 'Eiffel Tower' plot  (beta against psi) or a dotchart of beta if psi=FALSE.
##' @title Plot model parameters
##' @param x wordfish model parameters
##' @param pch plotting character: defaults to fat dots
##' @param psi whether to plot beta or beta against psi
##' @param ... Not used
##' @return Nothing. Used for plotting side-effect
##' @export
##' @method plot coef.wordfish
##' @author Will Lowe
plot.coef.wordfish <- function(x, pch=20, psi=TRUE, ...){  
  if (is.null(x$docs))
    stop(paste("Plotting word parameters in the multinomial parameterization",
               "probably won't\n  be very informative.  Try plotting the value",
               "of coef(model, form='poisson')"))
  
  ord <- order(x$words$beta)
  
  if (!psi){
    dotchart(x$words$beta[ord],
             labels=row.names(x$words)[ord],
             pch=pch, ...)
  } else {
    plot(x$words$beta[ord],
         x$words$psi[ord],
         type='n',
         xlab="Beta",
         ylab="Psi", ...)
    text(x$words$beta[ord],
         x$words$psi[ord],
         row.names(x$words)[ord])
  }
}

##' Summarise a wordfish model
##'
##' Summarises the fitted document positions theta (row scores) with confidence intervals
##' 
##' @title Summarise a wordfish model
##' @param object fitted wordfish model
##' @param level coverage as a proportion
##' @param ... Not used
##' @return A lightweight summary object containing the model call and position estimates
##' @export
##' @method summary wordfish
##' @author Will Lowe
summary.wordfish <- function(object, level=0.95, ...){
  m <- object
  z <- qnorm(1-(1-level)/2)
  pp <- data.frame(Estimate=m$theta,
                   "Std. Error"=m$se.theta,
                   Lower=(m$theta - m$se.theta*z),
                   Upper=(m$theta + m$se.theta*z))
  
  rownames(pp) <- m$docs
  ret <- list(model.call=m$call, scores=pp)
  class(ret) <- c('summary.wordfish', class(ret))
  return(ret)
}

##' Print a wordfish model summary
##'
##' Prints a wordfish model summary
##' @title Print a wordfish model summary
##' @param x wordfish model summary
##' @param digits how many significant digits to use
##' @param ... Not used
##' @return None. Used for printing side-effect
##' @export
##' @method print summary.wordfish
##' @author Will Lowe
print.summary.wordfish <- function(x, digits=max(3,getOption('digits')-3), ...){
  cat("Call:\n\t")
  print(x$model.call)
  cat("\nDocument Positions:\n")
  print(x$scores, digits=digits)
}

##' Predictions for new documents
##'
##' This function provides predicted values and uncertainty estimates for new documents
##' assuming that other model coefficients are known.  The multinomial formulation is
##' used to create analytic measures of uncertainty.  A better version of this function
##' would take into account uncertainty in the word parameter values when constructing
##' the latter.
##' 
##' @title Predictions for new documents. 
##' @param object fitted wordfish model 
##' @param newdata new data or the original data used to fit the model if NULL
##' @param se.fit whether standard errors are reported
##' @param interval whether confidence intervals are reported
##' @param level coverage as a proportion
##' @param ... Not used
##' @return a data.frame containing predicted positions for each document in newdata
##' @export
##' @method predict wordfish
##' @author Will Lowe
predict.wordfish <- function(object, newdata=NULL, se.fit=FALSE, interval=c("none", "confidence"), level=0.95, ...){
  m <- object
  if (is.null(newdata))
    newdata <- m$data
  interval <- match.arg(interval)
  if (interval != "none")
    se.fit=TRUE
  
  ## asymptotic standard errors
  mnll <- function(tt, b, p, y){
    linpred <- c(0, tt*b + p)
    dmultinom(y, size=sum(y), exp(linpred), log=TRUE)
  }
  V <- length(m$beta)
  new.beta <- (m$beta[2:V] - m$beta[1])
  new.psi <- (m$psi[2:V] - m$psi[1])
  
  preds <- rep(NA, nrow(newdata))
  preds.se <- rep(NA, nrow(newdata))
  for (i in 1:nrow(newdata)){
    preds[i] <- optimize(mnll, interval=c(-6,6), new.beta, new.psi, newdata[i,], maximum=TRUE)$maximum
    if (se.fit){
      neghess <- -hessian(mnll, preds[i], b=new.beta, p=new.psi, y=newdata[i,])
      invneghess <- solve(neghess)
      preds.se[i] <- sqrt(invneghess)
    }
  }
  if (se.fit){
    res <- data.frame(fit=preds, se.fit=preds.se)
    rownames(res) <- rownames(newdata)
    if (interval == "confidence"){
      z <- qnorm(1-(1-level)/2)
      res$lwr <- preds - z * preds.se
      res$upr <- preds + z * preds.se
    }
  } else {
    res <- preds
    names(res) <- rownames(newdata) ## we can at least label it
  }
  return(res)
}

##' Plot a wordfish model
##'
##' Plots document postions, associated measures of uncertainty, and optionally
##' true positions if they are available.
##' 
##' @title Plot a wordfish model 
##' @param x fitted wordfish model
##' @param truevals true positions for comparison
##' @param level coverage as a proportion
##' @param pch plotting character
##' @param ... parameters passed to dotchart
##' @return Nothing. Used for plotting side effect
##' @export
##' @method plot wordfish
##' @author Will Lowe
plot.wordfish <- function(x, truevals=NULL, level=0.95, pch=20, ...){
  ord      <- order(x$theta)
  theta    <- x$theta[ord]
  ci.theta <- (x$se.theta * qnorm(1-(1-level)/2))[ord]
  upper    <- theta + ci.theta
  lower    <- theta - ci.theta
  name.theta <- x$docs[ord]
  
  if (!is.null(truevals)){
    truevals <- truevals[ord]
    lims <- c(min(c(theta, truevals, lower)), max(c(theta, truevals, upper)))
    dotchart(theta, labels=name.theta, xlim=lims, pch=pch, ...)
    segments(lower, 1:length(theta), upper, 1:length(theta))
    points(truevals, 1:length(theta), col=rgb(139/255,0,0,0.75), pch=pch)
    title(paste('r =', format(cor(truevals, x$theta), digits=4)))
  } else {
    lims <- c(min(c(theta, lower)), max(c(theta, upper)))
    dotchart(theta, labels=name.theta, xlim=lims, pch=pch, ...)
    segments(lower, 1:length(theta), upper, 1:length(theta))
  }
}

##' Simulate wordfish data
##'
##' Simulates wordfish data where psi=2, alpha[1]=0, alpha[2:docs]=1,
##' theta consists of docs evenly spaced values with mean 0 and variance
##' 1, and beta consists vocab evenly spaced values between -1.5 and 1.5.
##' A simulated data matrix Y and these parameter values are returned a
##' named list.
##' 
##' @title simulate wordfish data 
##' @param docs how many documents to create
##' @param vocab how many words to create
##' @return a named list containing simulated data and parameter values
##' @export
##' @author Will Lowe
sim.wordfish <- function(docs=10, vocab=20){
  psi <- rep(2, vocab)
  alpha <- c(0, rep(1, docs-1))
  x <- 1:docs
  sdev <- sqrt(mean( (x - mean(x))**2))
  theta <- (x - mean(x)) / sdev
  beta <- seq(-1.5, 1.5, length.out=vocab)

  EY <- exp(outer(theta, beta) + outer(alpha, psi, "+"))
  Y <- matrix(rpois(rep(1, length(EY)), lambda=as.numeric(EY)),
              nrow=docs, ncol=vocab) ## columnwise out and in
  padder <- function(max.num, letter){
    sprintf(paste(letter, "%0", nchar(max.num), "d", sep=''), 1:max.num)
  }
  dimnames(Y) <- list(docs=padder(docs, "D"), words=padder(vocab, "W")) 
  list(Y=Y, psi=psi, beta=beta, alpha=alpha, theta=theta)
}

##' Fitted values for a wordfish model
##'
##' The function generates a matrix of expected counts the same size as the
##' data matrix used to fit the model.
##' 
##' @title Fitted values for a wordfish model 
##' @param object fitted wordfish model
##' @param ... Not used
##' @return matrix of expected counts
##' @export
##' @method fitted wordfish
##' @author Will Lowe
fitted.wordfish <- function(object, ...){
  m <- object
  n <- as.numeric(colSums(m$data))
  eta <- outer(m$theta, m$beta) + outer(m$alpha, m$psi, "+")
  mm <- exp(eta)
  ## mm <- apply(exp(eta), 1, function(x){ x / sum(x) })
  ## mm <- mm * outer(rep(1,length(m$beta)), n)
  colnames(mm) <- m$docs
  rownames(mm) <- m$words
  return(mm)
}

