##' Summarise a wordscores model
##'
##' Generates some minimal statistics about the document set from
##' which wordscores have been calculated.
##' 
##' @title Summarise a wordscores model 
##' @param object wordscores model
##' @param ... Not used
##' @return A data.frame with document statistics
##' @method summary wordscores
##' @export
##' @author Will Lowe
summary.wordscores <- function(object, ...){
  cat("Call:\n\t")
  print(object$call)
  
  cat("\nReference Document Statistics:\n\n")
  dd <- data.frame(Total=apply(object$data, 1, sum),
                   Min=apply(object$data, 1, min),
                   Max=apply(object$data, 1, max),
                   Mean=apply(object$data, 1, mean),
                   Median=apply(object$data, 1, median),
                   Score=object$theta)
  rownames(dd) <- rownames(object$data)
  print(dd, digits=3)
  invisible(dd)
}

##' Wordscores from a wordscores model
##'
##' Just lists the wordscores
##' @title Wordscores from a wordscores model
##' @param object a wordscores model
##' @param ... Not used
##' @return A set of scores for words
##' @export
##' @method coef wordscores
##' @author Will Lowe
coef.wordscores <- function(object, ...){
  return(object$pi)
}

##' Plot the wordscores from a model
##' 
##' Makes a dotchart of wordscores, order in size.
##' Not especially informative, except for viewing the document-driven
##' spikes described in Lowe 2008.
##' 
##' @title Plot wordscores 
##' @param x a wordscores model
##' @param ...  Graphics parameters passed on to dotchart
##' @return Nothing. Used for plotting side effect
##' @export
##' @method plot wordscores
##' @author Will Lowe
plot.wordscores <- function(x, ...){
  ord <- order(x$pi)
  dotchart(x$pi[ord], labels=names(x$pi)[ord], ...)
}

##' Score new documents with a wordscores model
##'
##' This function uses the wordscores from a fitted model
##' to estimate document positions for new documents, a.k.a.
##' virgin documents.  If newdata=NULL
##' the original 'reference' documents are used (their scores are not in
##' general recovered due to shrinkage and other factors...)
##'
##' Standard errors are computed as per the original paper, which may or may not
##' be sensible; the idea of standard errors in the absence of an explicit
##' probability model is unclear.
##' 
##' @title Score new document with a wordscores model 
##' @param object wordscores model
##' @param newdata A new word coutn matrix.  If NULL, the original matrix
##' @param rescale whether to apply the rescaling described in Laver Benoit and Garry 2003 which fixes the new documents' variance to be the same as the original documents. 
##' @param level Notional coverage as a proportion
##' @param ... Not used.
##' @return A data.frame containing predited new document positions
##' @export
##' @method predict wordscores
##' @author Will Lowe
predict.wordscores <- function(object, newdata=NULL, rescale=c('lbg', 'none'), level=0.95, ...){
  m <- object
  if (is.null(newdata))
    newd <- m$data
  else 
    newd <- newdata
  
  scorable <- which(colnames(newd) %in% names(m$pi))
  pi <- as.vector(m$pi[colnames(newd)[scorable]])

  cat(length(scorable), " of ", length(colnames(newd)), " words (",
      round(100*length(scorable)/length(colnames(newd)), 2), "%) are scorable\n\n",
      sep='')
  scorable.newd <- newd[,scorable,drop=FALSE]
  
  preds <- apply(scorable.newd, 1, function(x){ sum(x*pi)/sum(x) }) ## point estimate
  rowsum <- rowSums(scorable.newd) ## doc lengths
  preds.se <- rep(0, length(preds))
  for (i in 1:length(preds)){
    preds.se[i] <- sqrt(sum(scorable.newd[i,] * (preds[i] - pi)**2 / rowsum[i])) / sqrt(rowsum[i])
  }

  z <- qnorm(1-(1-level)/2)
  rs <- match.arg(rescale)
  if (rs=='lbg'){
    SDr <- sd(m$theta)
    Sv <- mean(preds)
    SDv <- ifelse(length(preds)<2, 0, sd(preds))
    mult <- ifelse(SDv==0, 0, SDr/SDv)
    re.theta <- (preds - Sv) * mult + Sv
    if (mult==0){
      ## corner case for no variance pointing out the bogosity of rescaling
      int.high <- preds + z * preds.se
      int.low <- preds - z * preds.se
    } else {
      int.high <- ((preds + z * preds.se) - Sv) * mult + Sv
      int.low <- ((preds - z * preds.se) - Sv) * mult + Sv
    }
    dd <- matrix(cbind(preds, preds.se, re.theta, int.low, int.high), ncol=5)
    colnames(dd) <- c("Score", "Std. Err.", "Rescaled", "Lower", "Upper")
  } else {
    dd <- matrix(cbind(preds, se.pres=preds.se), ncol=2)
    colnames(dd) <- c("Score", "Std. Err.")
  }
  rownames(dd) <- rownames(scorable.newd)
  
  print(dd, digits=3)
  invisible(as.data.frame(dd))
}

