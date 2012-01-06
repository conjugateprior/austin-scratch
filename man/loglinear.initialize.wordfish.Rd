\name{loglinear.initialize.wordfish}
\alias{loglinear.initialize.wordfish}
\title{Initialise wordfish parameters from a loglinear model}
\usage{
  loglinear.initialize.wordfish(wfm)
}
\arguments{
  \item{wfm}{word count matrix with documents as rows and
  words as columns}
}
\value{
  An initial set of wordfish model parameters
}
\description{
  Initialise wordfish parameters from a loglinear model
}
\details{
  This function uses a loglinear model as the basis for
  parameter estimates. First, we use loglin to fit a
  saturated model, adding a half to all counts to avoid
  numerical difficulties.  The log linear parameters
  consist of an intercept, the marginal effects 'docs' and
  'words', and a matrix of interaction terms 'docs.words'.
  The interaction matrix is then decomposed via SVD.
  Starting values for theta are taken from the first column
  of U and beta values from the first column of U (in the
  orientation given by R's svd function).  Starting values
  for alpha and psi are taken from the intercept and main
  effect parameters.  All parameters are then transformed
  so that theta is zero mean unit variance, alpha[1] is
  zero, and the intercept is absorbed into the main effects
  as per the usual wordfish parameterisation.
}
\author{
  Will Lowe
}

