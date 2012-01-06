\name{print.coef.wordfish}
\alias{print.coef.wordfish}
\title{Print model parameters}
\usage{
  \method{print}{coef.wordfish} (x,
    digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{wordfish model parameters from coef}

  \item{digits}{how many decinmal places to show}

  \item{...}{Not used}
}
\value{
  Nothing. Used for printing side-effect
}
\description{
  Print model parameters
}
\details{
  Prints out the model parameters in whatever
  parameterisation they were asked for in.
}
\author{
  Will Lowe
}

