\name{classic.wordfish}
\alias{classic.wordfish}
\title{Fit a wordfish model in S&P parameterisation}
\usage{
  classic.wordfish(wfm, dir, tol, sigma, params,
    verbose = FALSE)
}
\arguments{
  \item{wfm}{a word count matrix: documents are rows and
  columns are words}

  \item{dir}{identification constraint: theta[dir[1]] will
  be less than theta[dir[2]]}

  \item{tol}{the log likelihood amount less than which
  parameter estimation will stop}

  \item{sigma}{ridge regression regularization parameter}

  \item{params}{an initialized set of parameter to start
  estimation from}

  \item{verbose}{whether to give running commentary}
}
\value{
  a fitted wordfish model
}
\description{
  Fit a wordfish model S&P style
}
\details{
  This engine fits the Slapin and Proksch parameterisation
  of a regularized version of Goodman's RC model by
  coordinate ascent, alternating the fitting of row scores
  theta (with main effects alpha) and column scores beta
  (with main effects psi). This version of the model is
  identified by imposing 1) a zero mean and unit variance
  constraint on the row scores (theta), 2) a fixed zero for
  the first row effect (alpha), and 3) a directional
  constraint on row scores: theta[dir[1]] < theta[dir[2]].

  There are some minor differences to previous
  implemenentations: Unlike Goodman's RC model, the beta is
  the MAP estimated computed with a mean zero, standard
  deviation sigma Normally distrubuted prior. Unlike Slapin
  and Proksch's implementation, 'standard errors' for theta
  are computed from the Hessian, assuming other parameters
  fixed, rather than using a parametric bootstrap. This is
  quicker, asymptotically equivalent, and tends to give
  very similar results.
}
\author{
  Will Lowe
}

