\name{parametric.bootstrap.se}
\alias{parametric.bootstrap.se}
\title{Parametric bootstrap for a wordfish model}
\usage{
  parametric.bootstrap.se(m, L = 50,
    fit.fun = "classic.wordfish", verbose = FALSE)
}
\arguments{
  \item{m}{wordfish model}

  \item{L}{number of boostrap samples}

  \item{fit.fun}{engine used to fit the model.  Should be
  the same as the engine used to fit m}

  \item{verbose}{whether to generate a running commentary}
}
\value{
  For D documents, a D by L matrix of bootstrapped theta
  values (row scores)
}
\description{
  Parametric bootstrap for a wordfish model
}
\details{
  Runs a parametric bootstrap on a wordfish model, i.e.
  asuming that the model is true.  This is asymptotically
  equivalent to the Hessian-derived standard errors
  provided by wordfish.  The function creates L new word
  frequency matrices by sampling from the fitted model and
  returns a matrix of estimated thetas (row scores) To
  summarise these as standard errors take the standard
  deviation or compute quantiles of the returned matrix
  across rows.

  The default value of L is unrealistically small,
  particularly if quantile statistics are used to summarise
  the bootstrap samples.

  Note: With some probability words that occur in the
  original word frequency matrix will not appear in any
  bootstrapped word frequency matrix.  These are quietly
  discarded before fitting a model to the bootstrapped
  counts.
}
\author{
  Will Lowe
}

