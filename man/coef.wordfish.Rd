\name{coef.wordfish}
\alias{coef.wordfish}
\title{Show wordfish model parameters}
\usage{
  \method{coef}{wordfish} (object,
    form = c("poisson", "multinomial"), ...)
}
\arguments{
  \item{object}{fitted wordfish model}

  \item{form}{which parameterization to return}

  \item{...}{Not used.}
}
\value{
  a list containing two sets of word parameters and
  document main effects if form='poisson'
}
\description{
  Show wordfish model parameters
}
\details{
  Shows the parameters of a wordfish model in various
  equivalent parameterisations. Currently these are
  'poisson', the original Wordfish and RC model version,
  and 'multinomial', the equivalent model when row totals
  (document lengths) are assumed to be fixed.

  Note: This function does not return estimated document
  positions (row scores) theta. Use summary for that.
}
\author{
  Will Lowe
}

