\name{wordfish}
\alias{wordfish}
\title{Fit a wordfish model}
\usage{
  wordfish(wfm, dir = c(1, 10), tol = 1e-06, sigma = 3,
    params = NULL,
    init.fun = "loglinear.initialize.wordfish",
    fit.fun = "classic.wordfish", verbose = FALSE)
}
\arguments{
  \item{wfm}{a word count matrix, documents are rows and
  columns are words}

  \item{dir}{identification constraint: theta[dir[1]] will
  be less than theta[dir[2]]}

  \item{tol}{the log likelihood amount less than which
  parameter estimation will stop}

  \item{sigma}{is ridge regression regularization
  parameter}

  \item{params}{is an initialized set of parameter to start
  estimation from or NULL, which will trigger the
  application of init.fun}

  \item{init.fun}{function to generate an initial set of
  parameters. By default: loglinear.initialize.wordfish}

  \item{fit.fun}{the engine to use to fit the model. By
  default: classic.wordfish}

  \item{verbose}{whether to give running commentary on
  estimation}
}
\value{
  a fitted wordfish model
}
\description{
  Fit a wordfish model
}
\details{
  This engine fits the Slapin and Proksch parameterisation
  of a regularized version of Goodman's RC model by
  coordinate ascent, alternating the fitting of row scores
  theta (with main effects alpha) and column scores beta
  (with main effects psi). This version of the model is
  identified by imposing 1) zero mean and unit variance
  constrants on row scores (theta), 2) a fixed zero for the
  first row effect (alpha), and 3) a directional constraint
  on row scores: theta[dir[1]] < theta[dir[2]].

  There are some minor differences to previous
  implemenentations: Unlike Goodman's RC model, beta is
  ridge regularised with standard deviation sigma, there is
  no grand mean and no separate correlation parameter.
  Unlike Slapin and Proksch's implementation, 'standard
  errors' for theta are computed from the Hessian, assuming
  other parameters fixed rather than with a parametric
  bootstrap.
}
\author{
  Will Lowe
}

